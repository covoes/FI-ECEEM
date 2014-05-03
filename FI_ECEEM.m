function [bestPartition tFinal g gmmObjBest] = FI_ECEEM(data, constraints, configPrm)
%Feasible-Infeasible - Evolutionary Create & Eliminate EM algorithm
%
%%TODO refazer doc usando configPrm
%
%Parameters:
%data              : dataset to be clustered NxM where N is the number of objects and M the number of features
%constraints	     : pairwise constraints that should be respected Cx3 matrix with index to objects (starting from 1) in 1st and 2nd position
%					 and -1 for CL and 1 for ML in the 3rd position
%maxClusters       : maximum number of clusters to estimate
%sizePopulation    : number of individuals in the population
%maxEMIter         : maximum number of EM iterations for refinements in one EM execution
%fitnessFName      : name of fitness fuction to use. Possible are: 'mdl' (see fitnessFunc.m for details)
%maxKMSIter        : if the variance based splitting is used, define the number of kmeans iteration used for parameter estimation, also used in the initialization
%minSizePop		   : minimum size of feasible and infeasible populations.


if ischar(data) && strcmp(data,'debug')
	unittests();
	return
end




if ~isfield(configPrm,'DEBUG')
	configPrm.DEBUG=0;
end

if ~isfield(configPrm,'EXTRA_INFO')
	configPrm.EXTRA_INFO=0;
end
%DEBUG=fopen('debug.txt','w');

tIni = tic;

%%Control Variables
%tolerance to consider no improvement in EM
tolerance = 1e-5;

nonOptPrms = { 'maxClusters', 'sizePopulation',	'maxGenerations', 'maxGenWOImprov',
						   'maxEMIter', 'fitnessFName', 'maxKMSIter', 'minSizePop','regV'};
for f=nonOptPrms
	if ~isfield(configPrm,cell2mat(f))
		error('Field %s is missing from configPrm', cell2mat(f))
	end
end

configPrm.sizePopulationFeasible = configPrm.sizePopulation;
configPrm.sizePopulationInfeasible = configPrm.sizePopulation;

genWOImprov = 0;
lastBestFitness = 0;
conGraph = generate_constraint_graph(constraints, size(data,1));
[nChunklets,chunklets] = generate_chunklets(conGraph);
configPrm.minClusters = nChunklets;
staticSharedData = struct( 'constraints', constraints, 'conGraph', conGraph,...
	                         'data', data, 'nChunklets', nChunklets, 'chunklets', chunklets);

[Pfeas Pinfeas] = initialize_with_constraints(staticSharedData, configPrm);

g = 0;
while ~converged()
	g = g + 1;
	%TODO Trocar para parfor quando for rodar em paralelo
	newFeasibleSolutions = [];
	newInfeasibleSolutions = [];
	for i=1:length(Pfeas)
		indiv = Pfeas(i);
		if indiv.nClusters > 2
			indiv = clean_solution(indiv,staticSharedData, configPrm);
		end

		[indiv,gmmObj] = refinement(indiv, staticSharedData, configPrm);
		[feas infeas] = insert_individual_correct_pool(indiv, gmmObj, [], []);
		newFeasibleSolutions = [newFeasibleSolutions; feas];
		newInfeasibleSolutions = [newInfeasibleSolutions; infeas];

		[new_indiv,gmmObjNew] = feasible_mutation(indiv, gmmObj, staticSharedData, configPrm);
		[feas infeas] = insert_individual_correct_pool(new_indiv, gmmObjNew,[], []);
		newFeasibleSolutions = [newFeasibleSolutions; feas];
		newInfeasibleSolutions = [newInfeasibleSolutions; infeas];
	end

	if configPrm.DEBUG
		fprintf('%d - After FeasProc - %d feasPool %d infeasPool\n',g,length(newFeasibleSolutions),...
	          	length(newInfeasibleSolutions));
	end
	for i=1:length(Pinfeas)
		indiv = Pinfeas(i);
		[indiv, gmmObj] = refinement(indiv, staticSharedData, configPrm, 1);
		[new_indiv,gmmObj] = infeasible_mutation(indiv, gmmObj, staticSharedData, configPrm);
		[feas infeas] = insert_individual_correct_pool(new_indiv, gmmObj,[], []);
		newFeasibleSolutions = [newFeasibleSolutions; feas];
		newInfeasibleSolutions = [newInfeasibleSolutions; infeas];
	end

	if configPrm.DEBUG
		fprintf('%d - After InfeasProc - %d feasPool %d infeasPool\n',g,length(newFeasibleSolutions),...
	          	length(newInfeasibleSolutions));
	end
	feasiblePool = fitness_based_selection(Pfeas, newFeasibleSolutions, ...
		                                   configPrm.sizePopulationFeasible);

	infeasiblePool = constraint_based_selection(Pinfeas, newInfeasibleSolutions, ...
						                           configPrm.sizePopulationInfeasible);

	[Pfeas Pinfeas] = fill_pools_if_needed(staticSharedData, feasiblePool, ...
	                                     	infeasiblePool, configPrm);


	%fprintf('\n-%.2f (%.2f)\n',mean([ P(:).nClusters ]), std([P(:).nClusters]))
	if configPrm.EXTRA_INFO
		INFO_FIT(g,:) = [ Pfeas(:).fitness ];
		INFO_K(g,:) = [ Pfeas(:).nClusters ];
		INFO_FIT2(g,:) = [ Pinfeas(:).fitness ];
		INFO_K2(g,:) = [ Pinfeas(:).nClusters ];
		%save(EXTRA_INFO, 'INFO_FIT', 'INFO_K', 'INFO_FIT2', 'INFO_K2');
	end


end

tFinal = toc(tIni);

function converg = converged

	converg = 0;
	if isempty(Pfeas)
		return
	end

	[curBestFitness idx] = min([ Pfeas(:).fitness ]);
	if configPrm.DEBUG
		fprintf('%d -- MDL -- bst %.5f -- mean %.5f\n',g,curBestFitness, ...
			                             mean([Pfeas(:).fitness]));
		fprintf('%d -- VIO -- bst %.5f -- mean %.5f\n',g,min([Pinfeas(:).totPenalty]), ...
		                              	mean([Pinfeas(:).totPenalty]));
	end
	bestPartition = Pfeas(idx);
	bestPartition = clean_solution(bestPartition,staticSharedData, configPrm);
	[bestPartition, gmmObjBest] = refinement(bestPartition, staticSharedData, configPrm, 1);
	%test termination criterion
	if abs(lastBestFitness - curBestFitness) < tolerance
		genWOImprov = genWOImprov + 1;
		if genWOImprov == configPrm.maxGenWOImprov
			%no improvements in maxGenWOImprov iterations exit returning
			%the current best partition
			converg = 1;
		end
	else
		genWOImprov = 0;
		lastBestFitness = curBestFitness;
		converg = 0;
	end

	if g >= configPrm.maxGenerations
		converg = 1;
	end
end


end

function unittests
	testFIECEM_BestSolutionInfeasible;
	testFIECEM_BestSolutionFeasibleMultipleMappings;
	testFIECEM_BestSolutionFeasibleOneMapping;
end

function testFIECEM_BestSolutionFeasibleOneMapping
	data = mvnrnd([repmat([3 3],300,1); repmat([20 20], 300, 1)], [1 1]);
	constraints = [ 1 5 1; 5 90 1; 310 90 -1; 310 359 1];
	configPRM = struct('maxKMSIter',2,'maxClusters',5, 'sizePopulation',5, 'maxGenerations',5,...
		'maxGenWOImprov',2,'maxEMIter',3,'fitnessFName','mdl','minSizePop',2,'minClusters',2,...
		'maxInitTries',10, 'DEBUG',0, 'regV', 1e-2 );
	[bestPartition tFinal g] = FI_ECEEM(data, constraints, configPRM);
	[~,idx] = min(pdist2(bestPartition.mean,[3 3; 20 20]),[],2);
	mCorreta = [ 3 3; 20 20];
	cCorreta = [ 1 0 1; 1 0 1];
	assertElementsAlmostEqual(bestPartition.mean, mCorreta(idx,:), 'absolute',0.2)
	assertElementsAlmostEqual(bestPartition.covariance, cCorreta(idx,:), 'absolute',0.3)
	assertElementsAlmostEqual(bestPartition.mixCoef, [0.5 0.5], 'absolute',0.00001)
	assertTrue(isequal(bestPartition.classOfCluster(idx), [1;2]))
end


function testFIECEM_BestSolutionFeasibleMultipleMappings
	data = mvnrnd([repmat([3 3],300,1); ...
								 repmat([3 20], 300, 1); ...
	               repmat([20 20], 300, 1); ...
		             repmat([20 3], 300, 1)], [1 1]);
	constraints = [ 1 5 1; 5 90 1; 310 90 -1; 310 359 1; 610 810 1; 590 810 -1; ...
		              990 1200 1; 1000 500 1; 1000 5 -1; 1200 610 -1; 590 90 -1; ...
		              610 90 1; 310 1100 1; 500 310 1; 990 1000 1];
	configPRM = struct('maxKMSIter',2,'maxClusters',5, 'sizePopulation',5, 'maxGenerations',10,...
		'maxGenWOImprov',2,'maxEMIter',3,'fitnessFName','mdl','minSizePop',2,'minClusters',2,...
		'maxInitTries',10, 'DEBUG',0, 'regV', 1e-5 );
	[bestPartition  tFinal g] = FI_ECEEM(data, constraints, configPRM);
	info_individual(bestPartition)
	assertTrue(bestPartition.nClusters == 4,'Wrong number of clusters')
	mCorreta = [3 3; 3 20; 20 20; 20 3];
	cCorreta = [ 1 0 1; 1 0 1; 1 0 1; 1 0 1];
	cClCorreta = [1;2;1;2];
	[~,idx] = min(pdist2(bestPartition.mean,mCorreta),[],2);
	assertElementsAlmostEqual(bestPartition.mean, mCorreta(idx,:), 'absolute',0.2)
	assertElementsAlmostEqual(bestPartition.covariance, cCorreta(idx,:), 'absolute',0.5)
	assertElementsAlmostEqual(bestPartition.mixCoef, [0.25 0.25 0.25 0.25], 'absolute',0.00001)
	assertTrue(isequal(bestPartition.classOfCluster, cClCorreta(idx)))
end



function testFIECEM_BestSolutionInfeasible
	rng(42)
	data = [mvnrnd(repmat([3 3], 500,1), [0.01 0; 0 0.75]); ...
	        mvnrnd(repmat([5 3], 500,1), [0.02 0; 0 0.85])];
	data2 = mvnrnd(repmat([4 3], 500,1), [0.02 0; 0 0.85]);
	c1 = find(data(:,2) > 3);
	c2 = find(data(:,2) <= 3);
	data = [data;data2];
	c3 = 1001:1050;
	nConSamples = 20;
	conC1 = randsample(c1, nConSamples);
	conC2 = randsample(c2, nConSamples);
	conC3 = randsample(c3, nConSamples);
	constraints = zeros([nConSamples*3 3]);
	idx = 1;
	for i=1:(nConSamples-1)
		constraints((idx:idx+5),:) = [ conC1(i) conC1(i+1) 1; ...
			                             conC2(i) conC2(i+1) 1; ...
			                             conC3(i) conC3(i+1) 1; ...
                               	 	 conC1(i) conC2(i) -1; ...
                                   conC1(i) conC3(i) -1; ...
                                   conC2(i) conC3(i) -1];
		idx = idx + 6;
	end

	configPRM = struct('maxKMSIter',2,'maxClusters',30, 'sizePopulation',10, 'maxGenerations',15,...
		'maxGenWOImprov',4,'maxEMIter',4,'fitnessFName','mdl','minSizePop',2,'minClusters',5,...
		'maxInitTries',100, 'DEBUG',0, 'regV', 1e-5, 'EXTRA_INFO', 0 );
	[bestPartition  tFinal gi gmmObj] = FI_ECEEM(data, constraints, configPRM);
	%for k=1:bestPartition.nClusters
	%	sprintf('%d - %d\n',k, sum(gmmObj.clusterLabels==k))
	%end
	%info_individual(bestPartition)
	%plot_individual(bestPartition,data)
	assertTrue(length(unique(gmmObj.classLabels(conC1))) == 1)
	assertTrue(length(unique(gmmObj.classLabels(conC2))) == 1)
	assertTrue(length(unique(gmmObj.classLabels(conC3))) == 1)
	assertTrue(unique(gmmObj.classLabels(conC1)) ~= unique(gmmObj.classLabels(conC2)))
	assertTrue(unique(gmmObj.classLabels(conC1)) ~= unique(gmmObj.classLabels(conC3)))
	assertTrue(unique(gmmObj.classLabels(conC3)) ~= unique(gmmObj.classLabels(conC2)))

function plot_individual(indiv,data)
	addpath('../../outros/matlab_libs')
	figure;
	c1 = find(data(1:1000,2) > 3);
	c2 = find(data(1:1000,2) <= 3);
	c3 = 1001:1500;
	hold all;
	plot(data(c1,1), data(c1,2), '.b')
	plot(data(c2,1), data(c2,2), '.r')
	plot(data(c3,1), data(c3,2), '.g')
	covs = zeros([2 2 indiv.nClusters]);
	for k=1:indiv.nClusters
		covs(:,:,k) = squareformSymmetric(indiv.covariance(k,:));
	end
	plotGMM(indiv.mean(indiv.classOfCluster==1,:)', covs(:,:,indiv.classOfCluster==1), [ .8 0 0 ], 1);
	plotGMM(indiv.mean(indiv.classOfCluster==2,:)', covs(:,:,indiv.classOfCluster==2), [ 0 .8 0 ], 1);
	plotGMM(indiv.mean(indiv.classOfCluster==3,:)', covs(:,:,indiv.classOfCluster==3), [ 0 0 .8 ], 1);

	plot(data(c1,1), data(c1,2), '.b')
	plot(data(c2,1), data(c2,2), '.r')
	plot(data(c3,1), data(c3,2), '.g')

	plot(data(conC1,1), data(conC1,2), 'sb', 'MarkerSize', 12)
	plot(data(conC2,1), data(conC2,2), '*r', 'MarkerSize', 12)
	plot(data(conC3,1), data(conC3,2), 'hm', 'MarkerSize', 12)

end
end
