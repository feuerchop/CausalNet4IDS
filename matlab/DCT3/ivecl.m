function transformeddata=ivecl(data, type)
% PURPOSE:
% Transform a vector into a square symmetric matrix with ones on the main diagonal 
% INPUTS:
% data:  A kx1 vector to be transformed to asymmetric NxN matrix.  
%        k must ba solution to the equation k^1-k-2*k=0
% type:  String with values 'full' If you want a full matrix or 'lower' if
%        you want a lower triangular matrix. Default is full
% 
% OUTPUTS:
% transformeddata - a N by N symetric square matrix 

% it is a small modification of the ivech function of Kevin Sheppard

% Warning: vecl is suitable for correlations. If you are working with
% covariances consider using vech instead. For more info about vecl and
% vech operators consult:"Handbook of Matrices" by H Lutkepohl, ISBN
% 0-471-97015-8, John Wiley and sons

% Author: Manthos Vogiatzoglou, UoM, 2008

if nargin==1
    type='full';
end
if size(data,2)>1 && size(data,1)>1
    error('data is a column vector')
end
if size(data,1)==1 && size(data,2)>1
    display('The row vector imput is transformed to column vector')
    data=data';
end
k=size(data,1);
sizeout=(1+sqrt(1+8*k))/2;
transformeddata=zeros(sizeout);
index=1;

for i=1:sizeout
    for j=(i+1):sizeout
        if strcmp(type,'full')==1
        transformeddata(i,i)=1;
        transformeddata(j,i)=data(index);
        transformeddata(i,j)=data(index);
        transformeddata(sizeout,sizeout)=1;
        elseif strcmp(type,'lower')==1
        transformeddata(i,i)=1;
        transformeddata(j,i)=data(index);
        transformeddata(sizeout,sizeout)=1;
        end
        index=index+1;
    end
end

