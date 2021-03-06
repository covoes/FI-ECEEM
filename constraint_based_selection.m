function [infeasiblePool] = constraint_based_selection(P, Pmut,sizePop)
%CONSTRAINT_BASED_SELECTION Performs roulette selection filtering solutions that violates constraints
%
%Assumes that the intent is to minimize the fitness function
	fullPop = [ P(:); Pmut(:) ];
	if isempty(fullPop)
		infeasiblePool = [];
	else
		%perform mu-lambda selection
		fitness = [fullPop(:).totPenalty];
		fitness = max(fitness) - fitness + 1;
		probs = fitness ./ sum(fitness);
		survivors = roulette_without_reposition( probs, min(length(probs),sizePop) );
		infeasiblePool = fullPop(survivors);
	end
end
