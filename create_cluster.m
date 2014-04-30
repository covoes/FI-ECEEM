function	[indiv] = create_cluster(indiv, idxObj, data, chunklets, oldClusters)

	if ischar(indiv) && strcmp(indiv,'debug')
		unittests();
		return
	end
	variances = var(data);
	nClusters = indiv.nClusters;
	indiv.nClusters = nClusters + 1;
	indiv.mean(nClusters+1,:) = data(idxObj,:);
	indiv.covariance(nClusters+1,:) = squareformSymmetric( 0.1*diag(variances) );
	indiv.mixCoef(nClusters+1) = indiv.mixCoef(oldClusters(idxObj))/2;
	indiv.mixCoef(oldClusters(idxObj)) = indiv.mixCoef(oldClusters(idxObj))/2;
	indiv.classOfCluster(nClusters+1,:) = pickClass(idxObj, data, chunklets);
	indiv.fitness = NaN;
	indiv.totPenalty = NaN;
	assert (isequal(sort(unique(indiv.classOfCluster(:)))', 1:max(chunklets)))
end

function [closestClass] = pickClass(idxObj, data, chunklets)
	if chunklets(idxObj) > 0
		closestClass = chunklets(idxObj);
	else
		labeled = data(chunklets>0,:);
		classes = chunklets(chunklets>0);
		dists = pdist2(data(idxObj,:), labeled, 'euclidean');
		%select class of the closest labeled object
		[~,closestLabeled] = min(dists);
		closestClass = classes(closestLabeled);
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

