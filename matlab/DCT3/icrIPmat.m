function out=icrIPmat(data)
% This is the inverse of the crIPmat. Data is a square CxC matrix with
% zeros below the main antidiagonal.
% for example a 3x3 matrix with rows as follows
% data(1,:)=[2 3 4]
% data(2,:)=[5 6 0]
% data(3,:)=[7 0 0]
% is a type of matrix that this function is written for
% icrIPmat transforms this matrix to a vector row wise, by taking only the
% elements above and on the main antidiagonal
% thus for this matrix data, the output of the function is the following
% vector: out=[2 3 4 5 6 7]'. If data is a square cell array c, with c{i,j}
% a vector 1xn, icrIPmat creates n column vectors.
if iscell(data)==0
C=size(data,2);
out=ones(size(C+1));
for i=1:C+1
    out(i,i)=1;
    for j=i+1:C+1
    out(i,j)=data(i,j-i);
    out(j,i)=out(i,j);
    end
end
out=vecl(out);
else
    [K,L]=size(data{1,1}); M=max(K,L);
    C=size(data,2);
    out=zeros(.5*C*(C+1),M);
    dum=1;
    for i=1:C
        for j=1:C+1-i
            out(dum,:)=data{i,j};
            dum=dum+1;
        end
    end
end
            




