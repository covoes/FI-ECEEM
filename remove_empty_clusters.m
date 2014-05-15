function indiv = remove_empty_clusters(indiv, gmmObj)

	if ischar(indiv) && strcmp(indiv,'debug')
		unittests();
		return
	end

	oldClassCluster = unique(indiv.classOfCluster);
	idxUsed = unique(gmmObj.clusterLabels);
	idxRem = true([1 indiv.nClusters]);
	idxRem(idxUsed) = false;
	indiv = remove_clusters(find(idxRem), indiv);
	assert(isequal(sort(unique(indiv.classOfCluster)),sort(oldClassCluster)));
end

function unittests
	test_removeEmpty;
	test_doNotRemove;
end

function test_removeEmpty
	indiv = struct('nClusters', 10,'mean', rand(10,2), 'covariance', rand(10,3),...
	               'mixCoef', [1/10 0 0 0 6/10 0 0 0 3/10 0],...
		             'classOfCluster', [1 2 3 1 2 3 2 3 3 3] );
	gmmObj = struct('clusterLabels',[1 5 5 5 5 5 5 9 9 9]);
	novo = remove_empty_clusters(indiv, gmmObj);
	assertEqual(novo.nClusters, 3);
	assertEqual(novo.mean, indiv.mean([1 5 9],:));
	assertEqual(novo.covariance, indiv.covariance([1 5 9],:));
	assertElementsAlmostEqual(novo.mixCoef, indiv.mixCoef([1 5 9]));
	assertEqual(novo.classOfCluster, [1 2 3]);
end
	
function test_doNotRemove
	indiv = struct('nClusters', 5,'mean', rand(5,2), 'covariance', rand(5,3),...
	               'mixCoef', rand(5,1),'classOfCluster', [1 3 2 3 1] );
	indiv.mixCoef = indiv.mixCoef/sum(indiv.mixCoef);
	gmmObj = struct('clusterLabels',[ 1 2 1 2 3 3 2 1 3 4 5 5 4 4 5 5 4]);
	novo = remove_empty_clusters(indiv, gmmObj);
	assertEqual(novo, indiv);
end
