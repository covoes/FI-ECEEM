function [indiv] = feasible_mutation(indiv, gmmObj, sharedData, configPrm)
%FEASIBLE_MUTATION Performs mutation on the indivuals in the population

global DEBUG;

if ischar(indiv) && strcmp(indiv,'debug')
	unittests();
	return
end

data = sharedData.data
chunklets = sharedData.chunklets

nClusters = indiv.nClusters;
pElim = (nClusters-2)/(configPrm.maxClusters-2);
if (rand() > pElim)
	mutOpToApply = 'create';
else
	mutOpToApply = 'elim';
end

if DEBUG
	fprintf(DEBUG,'\n\n\n-----%s\n\n#MUTATION\nOLD INDIVIDUAL (%d):%s\n',mutOpToApply,i,info_individual(indiv));
end

posterior = gmmObj.posterior;
gauss = gmmObj.pdf;
oldIndiv = indiv;
if strcmp(mutOpToApply,'elim')
	%remove clusters
	lnTerms = bsxfun(@plus, log(eps+indiv.mixCoef(1:nClusters)), log(eps+gauss));
	%transform in positive so that higher values will have more probability to mutate
	valuesMutationClusters = -1*sum(posterior .* lnTerms, 1);
	[valuesSorted idxSorted] = sort( valuesMutationClusters );
	valuesMutationClusters( idxSorted ) = 1:nClusters;

	probs = valuesMutationClusters ./ sum(valuesMutationClusters);
	for c=1:length(unique(indiv.classOfCluster))
		if sum(indiv.classOfCluster == c) == 1
			probs(indiv.classOfCluster == c) = 0
		end
	end
	if all(probs==0)
		return
	end	
	z = randi([1 nClusters-2]);
	%chosen are the clusters selected for removal
	chosen = roulette_without_reposition( probs, z );
	survivors = setdiff(1:nClusters, chosen);
	indiv.mean = indiv.mean(survivors,:);
	indiv.covariance = indiv.covariance(survivors,:);
	indiv.mixCoef = indiv.mixCoef(survivors);
	%re-normalize
	indiv.mixCoef = indiv.mixCoef ./ nansum(indiv.mixCoef);
	indiv.nClusters = length(survivors);
else
	%create new clusters
	logsPost = log2(eps+posterior);
	objEntropies = -sum(posterior.*logsPost,2);

	[valuesSorted idxSorted] = sort( objEntropies );
	objEntropies( idxSorted ) = 1:length(objEntropies);

	probs = objEntropies./sum(objEntropies);
	z = randi([1 configPrm.maxClusters-nClusters]);
	%chosen are the objects that will be used to create clusters
	chosen = roulette_without_reposition( probs, z );
	variances = var(data);
	%recover the clusters that were most probable to generate these points
	oldClusters = gmmObj.clusterLabels;
	for nc=1:z
		indiv = create_cluster(indiv, chosen(nc), data, chunklets, oldClusters);
	end
end

if DEBUG
	assert(indiv.nClusters >= 2, 'Too few clusters')
	assert(indiv.nClusters <= configPrm.maxClusters, 'Too many clusters')
	fprintf(DEBUG,'#MUTATION\nNEW INDIVIDUAL (%d):%s\n',i,info_individual(indiv));
end


end


