function indiv = fix_class_unlabeled_clusters(indiv, gmmObj, sharedData)


	if ischar(indiv) && strcmp(indiv,'debug')
		unittests();
		return
	end


	labeled = sharedData.chunklets > 0;
	for k=1:indiv.nClusters
		objsInClusters = gmmObj.clusterLabels == k;
		if ~any(sharedData.chunklets(objsInClusters))
			indiv.classOfCluster(k) = closestClass(indiv.mean(k,:), ...
				                     sharedData.data(labeled,:), sharedData.chunklets(labeled));
		end
	end
end


function cls = closestClass(coord,labeled, classes)
	dists = pdist2(coord, labeled, 'euclidean');
	%select class of the closest labeled object
	[~,closestLabeled] = min(dists);
	cls = classes(closestLabeled);
end


function unittests
	test_fixClass;
end

function test_fixClass
	indiv = struct('nClusters',5, 'classOfCluster', [ 1 2 3 2 1], ...
		              'mean', [1 1; 6.9 6.9; 10.1 10.1; 2 2; 6 6 ]);
	gmmObj = struct('clusterLabels', [ 1 4 5 2 2 2 3 3 3 ]);
	sharedData = struct('chunklets', [ 1 0 0 2 2 2 3 3 3 ], ...
		                 'data', [1 1; 2 2; 6 6; 6 7; 7 6; 7 7; 6.5 7; 10 10; 11 10; 10 11]);
	new = fix_class_unlabeled_clusters(indiv, gmmObj, sharedData);
	assertTrue(isequal(new.classOfCluster, [1 2 3 1 2]));
end
