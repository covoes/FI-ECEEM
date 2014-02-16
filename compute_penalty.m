function [idxVio totPenalty] = compute_penalty(P, constraints, data)
%COMPUTE_PENALTY Given a population, compute a penalty score for each individual and return in 
%    idxVio a boolean vector of individuals with constraint violations

  if ischar(P) && strcmp(P,'debug')
		unittests()
		return
	end

	idxVio = logical(true([1 length(P)]));
	totPenalty = zeros([1 length(P)]);

	for i=1:length(P)
		totPenalty(i) = compute_penalty_model(P(i), data, constraints);
	end
  idxVio(totPenalty>0) = true;
end


function totPenalty = compute_penalty_model(individual, data, constraints, pdfs)
	if nargin < 4
		pdfs = computePosterior(individual, data);
	end
	totPenalty = 0;
	for c=1:size(constraints,1)
		idx1 = constraints(c,1);
		idx2  = constraints(c,2);
		typeCon = constraints(c,3);
		[~,k1] = max(pdfs(idx1,:));
		[~,k2] = max(pdfs(idx2,:));
		classK1 = individual.classOfCluster(k1);
		classK2 = individual.classOfCluster(k2);
		penalty = 0;
		if typeCon == 1 && classK1 ~= classK2
			%ML constraint being violated
			penalty = (1-pdfs(idx2, k1)) + (1-pdfs(idx1, k2));
		elseif typeCon == -1 && classK1 == classK2
			%CL constraind being violated
			penalty = pdfs(idx1, k1) +  pdfs(idx2, k2);
		end
		totPenalty = totPenalty + penalty;
	end
end



function unittests
	testComputePenaltyModel
end

function [data configPrm indiv con] = generateParamTesting
	data = [ 1 2; 1 3; 2 4; 5 5; 9 10; 10 11; 12 12];
	configPrm = struct('minClusters',2,'maxClusters',3, 'maxKMSIter',3, ...
		                 'sizePopFeasible',2, 'sizePopInfeasible', 2);
	indiv = struct('classOfCluster',[ 1 2 3]);
	con = [1 3 1; 6 7 1; 3 7 -1];
end


function testComputePenaltyModel

	[data configPrm indiv con] = generateParamTesting;

	testNoPenalty
	testMLPenalty
	testCLPenalty
	testMLPenalty_different_clusters_same_class
	testCLPenalty_different_clusters_same_class
	testMultiplePenalties

	function testNoPenalty
		gaussEasy = [ 1 0 0; 1 0 0; 1 0 0; 0 1 0; 0 0 1; 0 0 1; 0 0 1];
		out = compute_penalty_model(indiv, data, con, gaussEasy);
		assertTrue(out==0)
	end

	function testMLPenalty
		gaussEasy = [ 0.8 0.2 0; 1 0 0; 0.3 0.6 0.1; 0 1 0; 0 0 1; 0 0 1; 0 0 1];
		out = compute_penalty_model(indiv, data, con, gaussEasy);
		assertElementsAlmostEqual(out,1.5)
	end

	function testCLPenalty
		conN = [ 7 3 -1; 1 2 1; 6 7 -1];
		gaussEasy = [ 0.8 0.2 0; 1 0 0; 0.7 0.3 0; 0 1 0; 0 0 1; 0 0 1; 0.6 0 0.4];
		out = compute_penalty_model(indiv, data, conN, gaussEasy);
		assertElementsAlmostEqual(out,1.3)
	end

	function testCLPenalty_different_clusters_same_class
		conN = [ 4 5 -1 ];
		gaussEasy = [ 1 0 0; 1 0 0; 1 0 0; 0 1 0; 0 0 1; 0 0 1; 0 0 1];
		indiv = struct('classOfCluster',[ 1 2 2]);
		out = compute_penalty_model(indiv, data, conN, gaussEasy);
		assertElementsAlmostEqual(out,2)
	end

	function testMLPenalty_different_clusters_same_class
		gaussEasy = [ 0.8 0.2 0; 1 0 0; 0 0.9 0.1; 0 0.9 0.1; 0 0 1; 0 0 1; 0 0 1];
		indiv = struct('classOfCluster',[ 1 1 2]);
		out = compute_penalty_model(indiv, data, con, gaussEasy);
		assertElementsAlmostEqual(out,0)
	end

	function testMultiplePenalties
		indiv = struct('classOfCluster',[ 1 2 ]);
		conHard = [1 3 1; 6 7 1; 3 7 -1; 3 4 -1; 5 6 -1];
		gaussHard = [ 1 0; 1 0; 1 0; 0.55 0.45; 0 1; 0 1; 1 0];
		out = compute_penalty_model(indiv, data, conHard, gaussHard);
		assertElementsAlmostEqual(out, 7.55) ;
	end
end

