function [chosen] = roulette_without_reposition( probs, numToChoose )
%ROULETTE_WITHOUT_REPOSITION Performs roulette wheel selection without replacement

if ischar(probs) && strcmp(probs,'debug')
	unittests();
	return
end


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
end


%%TODO Terminar aqui
function test_roulette
	probs = [  0.1 0.3 0.6 ];
	freq = zeros([1 3]);
	nRpt = 10000;
	for r=1:nRpt
		counts = roulette_without_reposition(probs, 3);
		for c=1:3
			freq(c) = sum(counts==c)/nRpt;
		end
	end
	assertElementsAlmostEqual(probs, freq, 'absolute', 0.01)
end
