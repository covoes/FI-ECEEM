function [feasiblePool infeasiblePool] = constraint_based_selection(P, Pmut, constraints, feasiblePool, infeasiblePool, data)
%CONSTRAINT_BASED_SELECTION Performs (mu+lambda) selection filtering solutions that violates constraints
%
%Assumes that the intent is to minimize the fitness function
	global DEBUG;
	fullPop = [ P'; Pmut' ]; 
	%filter invalid solutions
	[idxVio,~] = check_constraints(fullPop, constraints, data);
	infeasiblePool = [infeasiblePool'; fullPop(idxVio)];
	fullPop(idxVio) = [];
	%perform mu-lambda selection
	fitness = [fullPop(:).fitness];
	[bestFit idxSorted] = sort(fitness,'ascend');
	P = fullPop(idxSorted(1:sizePopulation));
	feasiblePool = [feasiblePool'; P'];
end
