function [indiv,gmmObj] = infeasible_mutation(indiv, gmmObj, sharedData, configPrm)
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
penaltiesByCon = gmmObj.penalties;

nClusters = indiv.nClusters;
penalties = penaltiesByCon(:);
penalties = penalties./sum(penalties);
maxClusters = configPrm.maxClusters;
maxClustersToCreate = maxClusters - nClusters;
if maxClustersToCreate > 1
	constraintsToFix = randi(maxClustersToCreate,1) ;
else 
	constraintsToFix = 1;
end
chosen = roulette(penalties,constraintsToFix);

clusterLabels = gmmObj.clusterLabels;
for c=1:constraintsToFix
	con = constraints(chosen(c),:);
	po1 = clusterLabels(con(1));
	po2 = clusterLabels(con(2));
  classes = gmmObj.classLabels(con(1:2));

	if classes(1) == chunklets(con(1))
		assert(classes(2) ~= chunklets(con(2)))
		indiv = create_cluster(indiv, con(2), data, chunklets, clusterLabels) ;
	else
		assert(classes(1) ~= chunklets(con(1)))
		indiv = create_cluster(indiv, con(1), data, chunklets, clusterLabels) ;
	end
end

end
