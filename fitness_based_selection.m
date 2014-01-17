function [feasiblePool infeasiblePool] = fitness_based_selection(P, Pmut, constraints, feasiblePool, infeasiblePool)
%FITNESS_BASED_SELECTION Performs (mu+lambda) selection filtering solutions that violates constraints
%
%Assumes that the intent is to minimize the fitness function
	global DEBUG;
	fullPop = [ P'; Pmut' ]; 
	idxVio = check_constraints(fullPop, constraints, data);
	infeasiblePool = [infeasiblePool'; fullPop(idxVio)];
	fullPop(idxVio) = [];
	fitness = [fullPop(:).fitness];
	[bestFit idxSorted] = sort(fitness,'ascend');
	P = fullPop(idxSorted(1:sizePopulation));
	feasiblePool = [feasiblePool'; P'];
end
