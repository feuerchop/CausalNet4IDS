function out=crIPmat(data)
% PURPOSE:
% Transform a vector into a square matrix row-wise
% INPUTS:
% data:  A kx1 vector to be transformed to a  N-1xN-1 matrix.  
%        N must ba solution to the equation k^1-k-2*k=0

% OUTPUTS:
% out:  a N-1 by N-1 matrix. the elements above the main antidiagonal are 
%       filled by the elements of the vector data and all other elements are
%       equal to zeros. For example, if data is: data=[2 3 4 5 6 7]' then
%       the rows of out are given below
% out(:,1)=[2 3 4]
% out(:,2)=[5 6 0]
% out(:,3)=[7 0 0]

% Author: Manthos Vogiatzoglou, UoM, 2008
if iscell(data)==1
    out = data;
else
cols=size(data,2);
mat=ivecl(data(:,1));
[R, C]=size(mat);
if cols==1
    out=ones(C-1);
    for i=1:C-1
    out(:,i)=[diag(mat,i); zeros(i-1,1)];
    end
else
    out=cell(C-1);
    dum=0;
    for i=1:C-1
        for j=1:C-i
            dum=dum+1;
            out{i,j}=data(dum,:);
        end
    end
end
end
        
        
        

