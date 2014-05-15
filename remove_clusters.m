function 	indiv = remove_clusters(chosen, indiv)
%REMOVE_CLUSTERS Remove a list of clusters. 
% But if a cluster in the list is the only cluster of a class, it is not removed

	if ischar(chosen) && strcmp(chosen,'debug')
		unittests()
		return
	end

	chosen = filter_singleton_classes_clusters(chosen, indiv.classOfCluster);
	if isequal(chosen,[])
		return
	end
	oldIndiv = indiv;
	survivors = true([1 indiv.nClusters]);
	survivors(chosen) = false;
	indiv.nClusters = sum(survivors);
	indiv.mean = indiv.mean(survivors,:);
	indiv.covariance = indiv.covariance(survivors,:);
	indiv.mixCoef = indiv.mixCoef(survivors);
	indiv.classOfCluster = indiv.classOfCluster(survivors);
	%re-normalize
	indiv.mixCoef = indiv.mixCoef ./ nansum(indiv.mixCoef);

	assert(isequal(sort(unique(indiv.classOfCluster(:))), sort(unique(oldIndiv.classOfCluster(:)))))
end

function [filtered] = filter_singleton_classes_clusters(chosen, classOfCluster)
	filtered = chosen;
	classesOfChosen = classOfCluster(chosen);
	for c=1:max(classOfCluster)
		nClustersOfClass = sum(classOfCluster == c);
		chosenClustersOfClass = find(classesOfChosen == c);
		if nClustersOfClass == 1
			filtered(chosenClustersOfClass) = -1;
		elseif nClustersOfClass == length(chosenClustersOfClass)
			if nClustersOfClass == 0
				error('Class with no clusters')
			end
			filtered(chosenClustersOfClass) = -1;
			filtered(length(filtered)+1) = randsample(chosen(chosenClustersOfClass),1);
		end
	end
	filtered = filtered(filtered>-1);
end

function unittests
	test_filterSingletonClassesClusters
	test_removeClusters
end

function test_removeClusters
	indiv = struct('mean',[ 1 2; 3 4; 5 6; 7 8],'covariance',[ 1 0 1; 2 0 2; 3 0 3; 4 0 4;],...
		             'mixCoef', [ 0.25 0.25 0.25 0.25 ], 'classOfCluster', [ 1 1 2 2], ...
		             'nClusters', 4);
	nv = remove_clusters([2 3], indiv);
	assertElementsAlmostEqual(nv.mean, [1 2; 7 8])
	assertElementsAlmostEqual(nv.covariance, [1 0 1; 4 0 4])
	assertElementsAlmostEqual(nv.mixCoef, [0.5 0.5])
	assertElementsAlmostEqual(nv.classOfCluster, [1 2])
	assertElementsAlmostEqual(nv.nClusters, 2)
end

function test_filterSingletonClassesClusters
	chosen = [ 2 4 ];
	classOfClusters = [ 1 2 1 2];
	ret = filter_singleton_classes_clusters(chosen, classOfClusters);
	assertTrue(isequal(ret, 2) || isequal(ret,4));

	chosen = 2 ;
	classOfClusters = [ 1 2 1 ];
	ret = filter_singleton_classes_clusters(chosen, classOfClusters);
	assertTrue(isequal(ret, []))

	chosen = [ 1 4 ];
	classOfClusters = [ 1 2 1 3 4 3 ];
	ret = filter_singleton_classes_clusters(chosen, classOfClusters);
	assertTrue(isequal(ret, [1 4]))

	chosen = [ 1 3 4 6 ];
	classOfClusters = [ 1 2 1 3 4 3 ];
	ret = filter_singleton_classes_clusters(chosen, classOfClusters);
	assertTrue((ret(1) == 1 || ret(1) == 3) && (ret(2) == 4 || ret(2) == 6))
	assertTrue(length(ret) == 2)	
end


