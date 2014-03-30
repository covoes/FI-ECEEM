function [bestPartition EMSteps tFinal g] = FI_ECEEM(data, constraints, configPrm)
%Feasible-Infeasible - Evolutionary Create & Eliminate EM algorithm
%
%%TODO refazer doc usando configPRM
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


global EMSteps;

EXTRA_INFO=0;
%EXTRA_INFO='extraInfo.mat';


if ~isfield(configPrm,'DEBUG')
	configPrm.DEBUG=0;
end
%DEBUG=fopen('debug.txt','w');

tIni = tic;

%%Control Variables
%tolerance to consider no improvement in EM
tolerance = 1e-5;

nonOptPrms = {'maxClusters', 'sizePopulation', 'maxKMSIter', 'maxClusters', ...
	                'sizePopulation',	'maxGenerations', 'maxGenWOImprov', 'maxEMIter',...
                  'fitnessFName', 'maxKMSIter', 'minSizePop'};
for f=nonOptPrms
	if ~isfield(configPrm,cell2mat(f))
		error('Field %s is missing from configPrm', cell2mat(f))
	end
end

EMSteps = 0;
genWOImprov = 0;
lastBestFitness = 0;
conGraph = generate_constraint_graph(constraints, size(data,1));
[nChunklets,chunklets] = generate_chunklets(conGraph);
staticSharedData = struct( 'constraints', constraints, 'conGraph', conGraph,...
	                         'data', data, 'nChunklets', nChunklets, 'chunklets', chunklets);

[Pfeas Pinfeas] = initialize_with_constraints(staticSharedData, configPrm);

for g=1:configPrm.maxGenerations

	%TODO Trocar para parfor quando for rodar em paralelo
	newFeasibleSolutions = [];
	newInfeasibleSolutions = [];
	for i=1:length(Pfeas)
		indiv = Pfeas(i);

		[indiv,gmmObj] = refinement(indiv, staticSharedData, configPrm);
		[feas infeas] = insert_individual_correct_pool(indiv, gmmObj, staticSharedData,[], []);
		newFeasibleSolutions = [newFeasibleSolutions; feas];
		newInfeasibleSolutions = [newInfeasibleSolutions; infeas];

		new_indiv = feasible_mutation(indiv, gmmObj, staticSharedData, configPrm);
		[new_indiv,gmmObjNew] = refinement(new_indiv, staticSharedData, configPrm);
		[feas infeas] = insert_individual_correct_pool(new_indiv, gmmObjNew, staticSharedData,[], []);
		newFeasibleSolutions = [newFeasibleSolutions; feas];
		newInfeasibleSolutions = [newInfeasibleSolutions; infeas];
	end

	for i=1:length(Pinfeas)
		indiv = Pinfeas(i);
		[indiv, gmmObj] = refinement(indiv, staticSharedData, configPrm, 1);
		[new_indiv] = infeasible_mutation(indiv, gmmObj, staticSharedData, configPrm);
		[new_indiv, gmmObj] = refinement(new_indiv, staticSharedData, configPrm, 1);
		[feas infeas] = insert_individual_correct_pool(new_indiv, gmmObj, staticSharedData,[], []);
		newFeasibleSolutions = [newFeasibleSolutions; feas];
		newInfeasibleSolutions = [newInfeasibleSolutions; infeas];
	end

	feasiblePool = fitness_based_selection(Pfeas, newFeasibleSolutions);

	infeasiblePool = constraint_based_selection(Pinfeas, newInfeasibleSolutions);

	[feasiblePool infeasiblePool] = ...
			fill_pools_if_needed(feasiblePool, infeasiblePool, configPrm.minSizePop);

	Pfeas = feasiblePool(1:sizePopulation);
	Pinfeas = infeasiblePool(1:sizePopulation);

	if converged()
		return
	end


	%fprintf('\n-%.2f (%.2f)\n',mean([ P(:).nClusters ]), std([P(:).nClusters]))
	if EXTRA_INFO
		INFO_FIT(g,:) = [ Pfeas(:).fitness ];
		INFO_K(g,:) = [ Pfeas(:).nClusters ];
		INFO_FIT2(g,:) = [ Pinfeas(:).fitness ];
		INFO_K2(g,:) = [ Pinfeas(:).nClusters ];
		save(EXTRA_INFO, 'INFO_FIT', 'INFO_K', 'INFO_FIT2', 'INFO_K2');
	end

end

tFinal = toc(tIni);

function converged
	[curBestFitness idx] = min([ Pfeas(:).fitness ]);
	bestPartition = Pfeas(idx);

	%test termination criterion
	if abs(lastBestFitness - curBestFitness) < tolerance
		genWOImprov = genWOImprov + 1;
		if genWOImprov == maxGenWOImprov
			%no improvements in maxGenWOImprov iterations exit returning
			%the current best partition
			tFinal=toc(tIni);
			return ;
		end
	else
		genWOImprov = 0;
		lastBestFitness = curBestFitness;
		return ;
	end
end
end

function unittests
	testFIECEM;
end

function testFIECEM
	data = mvnrnd([repmat([3 3],100,1); repmat([20 20], 100, 1)], [1 1]);
	constraints = [ 1 2 1; 5 90 1; 110 90 -1; 122 159 1];
	configPRM = struct('maxKMSIter',2,'maxClusters',5, 'sizePopulation',5, 'maxGenerations',5,...
		'maxGenWOImprov',2,'maxEMIter',3,'fitnessFName','mdl','minSizePop',2,'minClusters',2,...
		'maxInitTries',10);
	[bestPartition EMSteps tFinal g] = FI_ECEEM(data, constraints, configPRM);
end
