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

posterior = gmmObj.posterior 
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
	indiv.nClusters = nClusters + z;
	%recover the clusters that were most probable to generate these points
	[~,oldClusters] = max(posterior(chosen,:), [], 2);
	for nc=1:z
		indiv.mean(nClusters+nc,:) = data(chosen(nc),:);
		indiv.covariance(nClusters+nc,:) = squareformSymmetric( 0.1*diag(variances) );
		indiv.mixCoef(nClusters+nc) = indiv.mixCoef(oldClusters(nc))/2;
		indiv.mixCoef(oldClusters(nc)) = indiv.mixCoef(oldClusters(nc))/2;
		indiv.classOfCluster(nClusters+nc,:) = pickClass(chosen(nc), data, chunklets)
	end
end

if DEBUG
	assert(indiv.nClusters >= 2, 'Too few clusters')
	assert(indiv.nClusters <= configPrm.maxClusters, 'Too many clusters')
	fprintf(DEBUG,'#MUTATION\nNEW INDIVIDUAL (%d):%s\n',i,info_individual(indiv));
end


end


function [closestClass] = pickClass(idxObj, data, chunklets)
	if chunklets(idxObj) > 0
		closestClass = chunklets(idxObj);
	else
		labeled = data(chunklets>0,:);
		classes = chunklets(chunklets>0);
		closestClass = NaN;
		dists = pdist2(data(idxObj,:), labeled, 'euclidean');
		dists = max(dists)-dists;
		dists = dists/sum(dists);
		chosen = roulette(dists, 1);
		closestClass = classes(chosen);
	end
end

function unittests
	testPickClass;
end

function testPickClass
	chunklets = [1 0 1 0 0 2 2];
	data = [ 1 2; 2 2; 1 2; 4 4; 80 80; 100 90; 120 100];
	assert(pickClass(1, data, chunklets) == 1)
	assert(pickClass(2, data, chunklets) == 1)
	assert(pickClass(3, data, chunklets) == 1)
	assert(pickClass(4, data, chunklets) == 1)
	assert(pickClass(5, data, chunklets) == 2)
	assert(pickClass(6, data, chunklets) == 2)
	assert(pickClass(7, data, chunklets) == 2)
end

