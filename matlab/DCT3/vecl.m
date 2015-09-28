function LTriangmat2vec=vecl(x)
% This function takes a square symmetric matrix of dimension NxN and stacks
% the elements below the main diagonal to a vector with size N(N-1)/2 x 1

% it is a minor modification of vech function of Kevin Sheppard, from ucsd
% garch toolbox

% Author: Manthos Vogiatzoglou, UoM, 2008

if size(x,1)~=size(x,2)
    display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    display('vecl has meaning only for square - symmetric matrices')
    display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    error('x should be square')
end
LTriangmat2vec=x(logical(tril(ones(size(x)),-1)));