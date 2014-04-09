function [nChunklets outChunklets] = generate_chunklets(conGraph)
%GENERATE_CHUNKLETS Obtain the chunklets basead on the constraints graph (only ML are needed).

  if ischar(conGraph) && strcmp(conGraph,'debug')
		unittests()
		return
	end
	MLs = conGraph>0;
	[nChunklets chunklets] = graphconncomp(MLs, 'Directed', false);
	validLabel = 1;
	outChunklets = chunklets;
	for i=1:nChunklets
		idx = find(chunklets == i);
		if numel(idx) == 1
			nChunklets = nChunklets - 1;
			outChunklets(idx) = 0;
		else
			outChunklets(idx) = validLabel;
			validLabel = validLabel + 1;
		end
	end
end

function unittests
	testGenerateChunklets
end

function testGenerateChunklets
	graph = generate_constraint_graph([1 3 1; 6 7 1; 4 5 -1], 7);
	[outNchunk chunklets] = generate_chunklets(graph);
	assertEqual(outNchunk, 2)
	assertEqual(chunklets, [1 0 1 0 0 2 2])

	graph = generate_constraint_graph([1 3 1; 6 7 1; 3 2 1], 7);
	[outNchunk chunklets] = generate_chunklets(graph);
	assertEqual(outNchunk, 2)
	assertEqual(chunklets, [1 1 1 0 0 2 2])
end

