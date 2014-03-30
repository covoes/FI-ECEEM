function [feasiblePool infeasiblePool] = fill_pools_if_needed(sharedData, ...
		feasiblePool, infeasiblePool,	configPrm)
%SORT_AND_FILL_POOLS_IF_NEEDED Makes sure that each pool has minSizePop individuals

infIndivToGenerate = configPrm.minSizePop - length(infeasiblePool);
feasIndivToGenerate = configPrm.minSizePop - length(feasiblePool); 

configPrm.sizePopulationFeasible = feasIndivToGenerate;
configPrm.sizePopulationInfeasible = infIndivToGenerate;
[feasiblePool2 infeasiblePool2] = initialize_with_constraints(sharedData, configPrm);

feasiblePool = [ feasiblePool(:); feasiblePool2(:) ];
infeasiblePool = [ infeasiblePool(:); infeasiblePool2(:) ];
