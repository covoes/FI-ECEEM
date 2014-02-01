function [bestPartition EMSteps tFinal g] = FI_ECEEM(data, constraints, maxClusters, sizePopulation, ...
	   	maxGenerations, maxGenWOImprov, maxEMIter, fitnessFName, maxKMSIter, minSizePop	)
%Feasible-Infeasible - Evolutionary Create & Eliminate EM algorithm
%
%
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


global EMSteps;
global DEBUG;

EXTRA_INFO=0;
%EXTRA_INFO='extraInfo.mat';


DEBUG=0;
%DEBUG=fopen('debug.txt','w');

tIni = tic;

%%Control Variables
%tolerance to consider no improvement in EM
tolerance = 1e-5;


EMSteps = 0;
genWOImprov = 0;
lastBestFitness = 0;

conGraph = generate_constraint_graph(constraints)

[Pfeas Pinfeas] = initialize_with_constraints(data, conGraph, maxClusters, sizePopulation, maxKMSIter);
for g=1:maxGenerations

	feasiblePool = []
	infeasiblePool = []
	Pfeas = refinement(Pfeas, data, maxEMIter, fitnessFName,maxKMSIter);
	[curBestFitness idx] = min([ Pfeas(:).fitness ]);

	bestPartition = P(idx);

	%test termination criterion
	if abs(lastBestFitness - curBestFitness) < tolerance
		genWOImprov = genWOImprov + 1;
		if genWOImprov == maxGenWOImprov
			%no improvements in maxGenWOImprov iterations exit returning
			%the current best partition
			tFinal=toc(tIni);
			g
			return
		end
	else
		genWOImprov = 0;
		lastBestFitness = curBestFitness;
	end

	Pmut = mutation(Pfeas, maxClusters, data);
	Pmut = refinement(Pmut, data, maxEMIter, fitnessFName, maxKMSIter);
	PmutInf = treatment(Pinfeas, data, maxEMIter);
	[feasiblePool infeasiblePool] = fitness_based_selection(Pfeas, Pmut, constraints, feasiblePool,
						infeasiblePool, data);
	[feasiblePool infeasiblePool] = constraint_based_selection(Pinfeas, PmutInf, constraints,
						feasiblePool, infeasiblePool, data);
	[feasiblePool infeasiblePool] = fill_pools_if_needed(feasiblePool, infeasiblePool, minSizePop);

	Pfeas = feasiblePool[1:sizePopulation];
	Pinfeas = infeasiblePool[1:sizePopulation];

	%fprintf('\n-%.2f (%.2f)\n',mean([ P(:).numClusters ]), std([P(:).numClusters]))
	if EXTRA_INFO
		INFO_FIT(g,:) = [ Pfeas(:).fitness ];
		INFO_K(g,:) = [ Pfeas(:).numClusters ];
		INFO_FIT2(g,:) = [ Pinfeas(:).fitness ];
		INFO_K2(g,:) = [ Pinfeas(:).numClusters ];
		save(EXTRA_INFO, 'INFO_FIT', 'INFO_K', 'INFO_FIT2', 'INFO_K2');
	end

end

tFinal = toc(tIni);

end
