function indiv = fix_class_unlabeled_clusters(indiv, gmmObj, sharedData)


	if ischar(indiv) && strcmp(indiv,'debug')
		unittests();
		return
	end

	oldIndiv = indiv;
	labeled = sharedData.chunklets > 0;
	nClustersPerClass = hist(indiv.classOfCluster, unique(indiv.classOfCluster));
	for k=1:indiv.nClusters
		objsInClusters = gmmObj.clusterLabels == k;
		if ~any(sharedData.chunklets(objsInClusters)) && ...
				nClustersPerClass(indiv.classOfCluster(k)) > 1
			indiv.classOfCluster(k) = closestClass(indiv.mean(k,:), ...
				                     sharedData.data(labeled,:), sharedData.chunklets(labeled));
		end
	end
	assert(isequal(sort(unique(indiv.classOfCluster(:))), ...
		             sort(unique(oldIndiv.classOfCluster(:)))))
end


function cls = closestClass(coord,labeled, classes)
	dists = pdist2(coord, labeled, 'euclidean');
	%select class of the closest labeled object
	[~,closestLabeled] = min(dists);
	cls = classes(closestLabeled);
end


function unittests
	test_fixClass;
	test_doNotRemoveSingleton;
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

function test_doNotRemoveSingleton
	%caso de soluções infeasible em que o objeto rotulado para um chunklet pode estar
	%em um grupo de outra classe, deve-se tomar cuidado adicional para não remover o 
	%único cluster daquela classe
	indiv = struct('nClusters',5, 'classOfCluster', [ 1 2 3 2 1], ...
		              'mean', [1 1; 6.9 6.9; 10.1 10.1; 2 2; 6 6 ]);
	gmmObj = struct('clusterLabels', [ 1 4 5 2 2 2 3 3 3 ]);
	sharedData = struct('chunklets', [ 1 0 3 2 2 2 0 0 0 ], ...
		                 'data', [1 1; 2 2; 6 6; 6 7; 7 6; 7 7; 6.5 7; 10 10; 11 10; 10 11]);
	new = fix_class_unlabeled_clusters(indiv, gmmObj, sharedData);
	assertEqual(new.classOfCluster, [1 2 3 1 1]);
end
