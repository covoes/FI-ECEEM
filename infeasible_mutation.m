function [indiv,gmmObj] = infeasible_mutation(indiv, gmmObj, sharedData, configPrm)
%FEASIBLE_MUTATION Performs mutation on the indivuals in the population


data = sharedData.data;
chunklets = sharedData.chunklets;

penaltiesByObj = gmmObj.penalties;

nClusters = indiv.nClusters;
penalties = penaltiesByObj(:);
penalties = penalties./sum(penalties);
maxClusters = configPrm.maxClusters;
maxClustersToCreate = maxClusters - nClusters;
if maxClustersToCreate > 1
	objsToFix = randi(maxClustersToCreate,1) ;
elseif maxClustersToCreate == 1
	objsToFix = 1;
else
	%maximum number of clusters reached
	%chosen are the clusters selected for removal
	chosen = randi(nClusters, [1 nClusters-configPrm.minClusters]);
	indiv = remove_clusters(chosen, indiv);
	return
end
chosen = roulette_without_reposition(penalties,objsToFix);

clusterLabels = gmmObj.clusterLabels;
for c=1:objsToFix
	indiv = create_cluster(indiv, chosen(c), data, chunklets, clusterLabels) ;
end

if indiv.nClusters > 2
	indiv = clean_solution(indiv,sharedData, configPrm);
end
[indiv, gmmObj] = refinement(indiv, sharedData, configPrm, 1);

end
