function [Pfeas Pinfeas] = insert_individual_correct_pool(individual, gmmObj, ...
		                                                      Pfeas, Pinfeas)
%INSERT_INDIVIDUAL_CORRECT_POOL Inserts a new individual in the feasible or infeasible pool
	if gmmObj.isFeasible
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
