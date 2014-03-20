function [feasiblePool infeasiblePool] = fill_pools_if_needed(feasiblePool, infeasiblePool,
						minSizePop)
%SORT_AND_FILL_POOLS_IF_NEEDED Makes sure that each pool has minSizePop individuals
%
%If the feasiblePool has fewer than minSizePop individuals, new individuals are generated
% using repair mechanisms in the individuals on the infeasiblePool.
%For the infeasiblePool new individuals are generated randomly.

infIndivToGenerate = minSizePop - len(infeasiblePool)

if infIndivToGenerate > 0
	newIndiv = initialize(data, maxClusters, infIndivToGenerate, maxKMSIter)
end
needFeasible = len(feasiblePool) < minSizePop
