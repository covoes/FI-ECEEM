function new = clean_solution(indiv, sharedData,configPrm)
%CLEAN_SOLUTION Remove empty clusters and fix classes for clusters with only unlabeleled objects
		[new, gmmObj] = refinement(indiv, sharedData, configPrm, 1);
		new = fix_class_unlabeled_clusters(new, gmmObj, sharedData);
		new = remove_empty_clusters(new, gmmObj);
		new = refinement(new, sharedData, configPrm, 1);
end

