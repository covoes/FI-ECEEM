function [infeasiblePool] = constraint_based_selection(P, Pmut)
%CONSTRAINT_BASED_SELECTION Performs (mu+lambda) selection filtering solutions that violates constraints
%
%Assumes that the intent is to minimize the fitness function
	fullPop = [ P(:); Pmut(:) ];
	%perform mu-lambda selection
	fitness = [fullPop(:).totPenalty];
	[~,idxSorted] = sort(fitness,'ascend');
	infeasiblePool = fullPop(idxSorted);
end
