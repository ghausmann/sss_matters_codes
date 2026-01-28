function [ A,ind ] = takecols( ten,cols )
%extract cols from a tensor class ten with symmetric columns.
%improve by mex
%
% © Copyright, Oren Levintal, June 13, 2016.
% bug corrected on Jan 11, 2018 for the cases of one or zero nonzero values in ten

m=length(cols);
temp=zeros(ten.tsize(2),1);
temp(cols)=1:m;
[i,j,vals]=tfind(ten);
j=temp(j);

n_vals=numel(i);
if n_vals==1
    j=j(:)';
end

if n_vals>0
    prodj=prod(j,2);
else
    prodj=0;
end

if nnz(prodj)>0
    A=sptensor(i(prodj~=0),j(prodj~=0,:),vals(:,prodj~=0),ten.tsize(1),repmat(m,1,length(ten.tsize)-1));
else
    A=sptensor(ten.tsize(1),repmat(m,1,length(ten.tsize)-1),size(ten.vals,1));
end

ind=find(prodj);

end

