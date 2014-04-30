function [feasiblePool] = fitness_based_selection(P, Pmut, sizePop)
%FITNESS_BASED_SELECTION Performs (mu+lambda) selection filtering solutions that violates constraints
%
%Assumes that the intent is to minimize the fitness function
	fullPop = [ P(:); Pmut(:) ]; 
	if isempty(fullPop)
		feasiblePool = [];
	else
		fitness = [fullPop(:).fitness];
		[~, idxSorted] = sort(fitness,'ascend');
		feasiblePool = fullPop(idxSorted(1:min(length(idxSorted),sizePop)));
	end
end
