function indiv = remove_empty_clusters(indiv, gmmObj)

	if ischar(indiv) && strcmp(indiv,'debug')
		unittests();
		return
	end

	oldClassCluster = unique(indiv.classOfCluster);
	idxUsed = unique(gmmObj.clusterLabels);
	indiv.mean = indiv.mean(idxUsed,:);
	indiv.covariance = indiv.covariance(idxUsed,:);
	indiv.mixCoef = indiv.mixCoef(idxUsed);
	indiv.classOfCluster = indiv.classOfCluster(idxUsed);
	indiv.nClusters = length(idxUsed);
	assert(isequal(unique(indiv.classOfCluster),oldClassCluster));
end

function unittests
	test_removeEmpty;
	test_doNotRemove;
end

function test_removeEmpty
	indiv = struct('nClusters', 10,'mean', rand(10,2), 'covariance', rand(10,3),...
	               'mixCoef', rand(10,1),'classOfCluster', [1 2 2 2 2 2 2 3 3 3] );
	gmmObj = struct('clusterLabels',[ 1 3 10]);
	novo = remove_empty_clusters(indiv, gmmObj);
	assertTrue(novo.nClusters == 3);
	assertTrue(isequal(novo.mean, indiv.mean([1 3 10],:)));
	assertTrue(isequal(novo.covariance, indiv.covariance([1 3 10],:)));
	assertTrue(isequal(novo.mixCoef, indiv.mixCoef([1 3 10])));
	assertTrue(isequal(novo.classOfCluster, [1 2 3]));
end
	
function test_doNotRemove
	indiv = struct('nClusters', 5,'mean', rand(5,2), 'covariance', rand(5,3),...
	               'mixCoef', rand(5,1),'classOfCluster', [1 3 2 3 1] );
	gmmObj = struct('clusterLabels',[ 1 2 1 2 3 3 2 1 3 4 5 5 4 4 5 5 4]);
	novo = remove_empty_clusters(indiv, gmmObj);
	assertTrue(isequal(novo, indiv));
end
