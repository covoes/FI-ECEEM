function run_tests(funcao)
addpath('~/workspace/codigos/outros/matlab_libs/matlab_xunit_3_1_1/xunit/')
addpath('~/workspace/codigos/proprios/utils/')
warning off stats:gmdistribution:FailedToConverge
warning off stats:kmeans:FailedToConverge
warning off stats:kmeans:EmptyCluster
warning on verbose
fh = str2func(funcao);
fh('debug')
