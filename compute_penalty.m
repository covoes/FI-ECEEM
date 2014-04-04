function [isInfeasible,totPenalty,penaltyByObj] = compute_penalty(individual, chunklets, ...
	                                   	post)
%COMPUTE_PENALTY Compute a penalty score for a individual 

	if ischar(individual) && strcmp(individual,'debug')
		unittests()
		return
	end
	labeled = find(chunklets>0);
	penaltyByObj = zeros([1 length(chunklets)]) ;
	for l=labeled
		[~,k] = max(post(l,:));
		classK = individual.classOfCluster(k);
		penalty = 0;
		if classK ~= chunklets(l)
			validClusters = individual.classOfCluster == chunklets(l);
			penalty = 1-max(post(l, validClusters));
		end
		penaltyByObj(l) = penalty;
	end
	totPenalty = sum(penaltyByObj);
	isInfeasible = totPenalty > 0;
end

function unittests
	testComputePenaltyModel
end

function [chunklets configPrm indiv] = generateParamTesting
	data = [ 1 2; 1 3; 2 4; 5 5; 9 10; 10 11; 12 12];
	configPrm = struct('minClusters',2,'maxClusters',3, 'maxKMSIter',3, ...
		                 'sizePopFeasible',2, 'sizePopInfeasible', 2);
	indiv = struct('classOfCluster',[ 1 2 3]);
	chunklets = [ 1 0 1 0 0 3 3 ];
end


function testComputePenaltyModel

	[chunklets configPrm indiv ] = generateParamTesting;

	testNoPenalty
	testPenalty
	testNoPenalty_different_clusters_same_class
	testPenalty_different_clusters_and_class
	testMultiplePenalties

	function testNoPenalty
		gaussEasy = [ 1 0 0; 1 0 0; 1 0 0; 0 1 0; 0 0 1; 0 0 1; 0 0 1];
		[isInf,tot,byobj]  = compute_penalty(indiv, chunklets,  gaussEasy);
		assertTrue(~isInf)
		assertTrue(tot==0)
		assert(all(byobj==0))
	end

	function testPenalty
		gaussEasy = [ 0.8 0.2 0; 1 0 0; 0.3 0.6 0.1; 0 1 0; 0 0 1; 0 0 1; 0 0 1];
		[isInf,tot,byobj] = compute_penalty(indiv, chunklets,  gaussEasy);
		assertTrue(isInf)
		assertTrue(byobj(3) == 0.7)
		assertTrue(all(byobj([1:2 4:7]) == 0))
		assertElementsAlmostEqual(tot,0.7)
	end

	function testNoPenalty_different_clusters_same_class
		chunklets = [1 0 1 2 2 2 2];
		gaussEasy = [ 1 0 0; 1 0 0; 1 0 0; 0 1 0; 0 0 1; 0 0 1; 0 0 1];
		indiv = struct('classOfCluster',[ 1 2 2]);
		[isInf,out] = compute_penalty(indiv, chunklets,  gaussEasy);
		assertElementsAlmostEqual(out,0)
		assertTrue(~isInf)
	end

	function testPenalty_different_clusters_and_class
		chunklets = [ 1 1 1 1 2 2 2 ];
		gaussEasy = [ 0.8 0.2 0; 1 0 0; 0 0.9 0.1; 0 0.1 0.9; 0 0 1; 0 0 1; 0 0 1];
		indiv = struct('classOfCluster',[ 1 1 2]);
		[isInf,tot,byobj] = compute_penalty(indiv, chunklets,  gaussEasy);
		assertElementsAlmostEqual(tot,0.9)
		assertElementsAlmostEqual(byobj, [0 0 0 0.9 0 0 0])
		assertTrue(isInf)
	end

	function testMultiplePenalties
		indiv = struct('classOfCluster',[ 1 2 ]);
		chunklets = [ 1 0 1 1 0 1 1 ];
		gaussHard = [ 1 0; 1 0; 1 0; 0.45 0.55; 0 1; 0 1; 1 0];
		[isInf,tot,byobj] = compute_penalty(indiv, chunklets,  gaussHard);
		assertElementsAlmostEqual(tot, 1.55) ;
		assertElementsAlmostEqual(byobj,[0,0,0,0.55,0,1,0])
		assertTrue(isInf)
	end
end

