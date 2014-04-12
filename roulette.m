function [chosen] = roulette( probs, numToChoose )
%ROULETTE Performs roulette wheel selection with replacement


if ischar(probs) && strcmp(probs,'debug')
	unittests();
	return
end

chosen = zeros(numToChoose,1);	

for i = 1:numToChoose
	draw = rand;
	acc = probs(1);
	chosen(i) = 1;
	while (draw > acc)
		chosen(i) = chosen(i) + 1;
		acc = acc + probs(chosen(i)); 
	end
end

end

function unittests
	test_roulette;
end

function test_roulette
	probs = [  0.1 0.3 0.6 ];
	freq = zeros([1 3]);
	nRpt = 10000;
	counts = roulette(probs, nRpt);
	for c=1:3
		freq(c) = sum(counts==c)/nRpt;
	end
	assertElementsAlmostEqual(probs, freq, 'absolute', 0.01)
end
