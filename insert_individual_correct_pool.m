function [Pfeas Pinfeas] = insert_individual_correct_pool(individual, gmmObj, sharedData, ...
		                                                      Pfeas, Pinfeas)
%INSERT_INDIVIDUAL_CORRECT_POOL Inserts a new individual in the feasible or infeasible pool
	if isFeasible(individual, sharedData.conGraph, sharedData.data, gmmObj)
		idx = length(Pfeas);
		if idx == 0
			Pfeas = individual;
		else
			Pfeas(idx+1) = individual;
		end
	else
		idx =  length(Pinfeas);
		if idx == 0
			Pinfeas = individual;
		else
			Pinfeas(idx+1) = individual;
		end
	end
end

function isF = isFeasible(individual, conGraph, data,gmmObj)
	[i j s] = find(conGraph);
	constraints = [ i j s ];
	idxVio = compute_penalty(individual, data, constraints, posterior(gmmObj, data));
	isF = idxVio == false;
end


