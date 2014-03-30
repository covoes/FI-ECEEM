function [feasiblePool] = fitness_based_selection(P, Pmut)
%FITNESS_BASED_SELECTION Performs (mu+lambda) selection filtering solutions that violates constraints
%
%Assumes that the intent is to minimize the fitness function
	fullPop = [ P(:); Pmut(:) ]; 
	if length(fullPop) == 0
		feasiblePool = [];
	else
		fitness = [fullPop(:).fitness];
		[bestFit idxSorted] = sort(fitness,'ascend');
		feasiblePool = fullPop(idxSorted);
	end
end
