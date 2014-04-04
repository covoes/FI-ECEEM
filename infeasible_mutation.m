function [indiv] = infeasible_mutation(indiv, gmmObj, sharedData, configPrm)
%FEASIBLE_MUTATION Performs mutation on the indivuals in the population

DEBUG = configPrm.DEBUG;

if ischar(indiv) && strcmp(indiv,'debug')
	unittests();
	return
end

constraints = sharedData.constraints;
data = sharedData.data;
nChunklets = sharedData.nChunklets;
chunklets = sharedData.chunklets;
conGraph = sharedData.conGraph;

pdfs = gmmObj.posterior;
totPenalty = indiv.totPenalty;
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
	return
end
chosen = roulette(penalties,objsToFix);

clusterLabels = gmmObj.clusterLabels;
for c=1:objsToFix
	indiv = create_cluster(indiv, chosen(c), data, chunklets, clusterLabels) ;
end

end
