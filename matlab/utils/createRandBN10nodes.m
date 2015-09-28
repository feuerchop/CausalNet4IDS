function bnet = createRandBN10nodes(base)
% setup nodes
A=1; B=2; C=3; D=4; E=5; F=6; G=7; H=8; I=9; J=10;
n = 10;

% setup graph 
dag = zeros(n);
dag(A, [B, C, D]) = 1;
dag(B, E) = 1;
dag(C, F) = 1;
dag(D, F) = 1;
dag(E, [G, H]) = 1;
dag(F, [H, J]) = 1;
dag(G, I) = 1;
dag(H, I) = 1;

% node sizes - all cts nodes are scalar, all discrete nodes are binary
ns = ones(1, n);
bnet = mk_bnet( dag, ns, 'names', {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J'}, 'discrete', [] );
bnet.CPD{A} = gaussian_CPD( bnet, A, 'mean', base*(randn(1)) , 'cov', base*(randn(1)));
bnet.CPD{B} = gaussian_CPD( bnet, B, 'mean', base*(randn(1)), 'cov', base*(randn(1)), 'weights', base*(randn(1)));
bnet.CPD{C} = gaussian_CPD( bnet, C, 'mean', base*(randn(1)), 'cov', base*(randn(1)), 'weights', base*(randn(1)));
bnet.CPD{D} = gaussian_CPD( bnet, D, 'mean', base*(randn(1)), 'cov', base*(randn(1)), 'weights', base*(randn(1)));
bnet.CPD{E} = gaussian_CPD( bnet, E, 'mean', base*(randn(1)), 'cov', base*(randn(1)), 'weights', base*(randn(1)));
bnet.CPD{F} = gaussian_CPD( bnet, F, 'mean', base*(randn(1)), 'cov', base*(randn(1)), 'weights', [base*(randn(1)), base*(randn(1))]);
bnet.CPD{G} = gaussian_CPD( bnet, G, 'mean', base*(randn(1)), 'cov', base*(randn(1)), 'weights', base*(randn(1)));
bnet.CPD{H} = gaussian_CPD( bnet, H, 'mean', base*(randn(1)), 'cov', base*(randn(1)), 'weights', [base*(randn(1)), base*(randn(1))]);
bnet.CPD{I} = gaussian_CPD( bnet, I, 'mean', base*(randn(1)), 'cov', base*(randn(1)), 'weights', [base*(randn(1)), base*(randn(1))]);
bnet.CPD{J} = gaussian_CPD( bnet, J, 'mean', base*(randn(1)), 'cov', base*(randn(1)), 'weights', base*(randn(1)));
end