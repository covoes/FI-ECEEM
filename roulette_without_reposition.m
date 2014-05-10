function [chosen] = roulette_without_reposition( probs, numToChoose )
%ROULETTE_WITHOUT_REPOSITION Performs roulette wheel selection without replacement

if ischar(probs) && strcmp(probs,'debug')
	unittests();
	return
end

assert(numToChoose <= sum(probs>0))

chosen = zeros(numToChoose,1);

for i = 1:numToChoose
	selected = roulette( probs, 1 );
	probs = probs / (1-probs(selected));
	probs( selected ) = 0;
	chosen(i) = selected;
end

end

function unittests
	test_roulette;
	test_rouletteWithZero;
end


function test_roulette
	probs = [  0.1 0.3 0.6 ];
	assertElementsAlmostEqual(probs, roda_roleta(probs), 'absolute', 0.01)
end

function test_rouletteWithZero
	probs = [ 0.1 0.4 0 0.02 0 0.01 0 0.18 0 0 0.09 0.2 ];
	assertElementsAlmostEqual(probs, roda_roleta(probs), 'absolute', 0.01)
end

function [freq] = roda_roleta(probs)
	freq = zeros([1 length(probs)]);
	nRpt = 10000;
	for r=1:nRpt
		counts = roulette_without_reposition(probs, sum(probs>0));
		for c=1:length(probs)
			freq(c) = freq(c) + (counts(1)==c);
		end
	end
	freq = freq/nRpt;
end
