function [ A,ind ] = takerows( ten,rows )
%extract rows from tensor ten.
%
% © Copyright, Oren Levintal, June 13, 2016.
% bug corrected on Jan 1, 2018 for the case of zero nonzero values in A
% bug corrected on Jan 29, 2018 for the case of an empty tensor

if ten.tsize(1)==0 && ~isempty(rows)
    error('trying to extract rows from an empty tensor')
elseif ten.tsize(1)==0 && isempty(rows)
    if nargout==1
        A=ten;
    else
        A=ten;
        ind=[];
    end
else
    if nargout==1 %mex version
        A=sptensor;
        IA=ten.ptr;
        JA=ten.cols;
        NA=ten.vals;
        l=ten.tsize(1);
        takel=intarray(rows);
        [IC,JC,NC] = takerows_mex(IA,JA,NA,l,takel);
        newl=length(takel);
        A.ptr=IC(1:newl+1);
        A.cols=JC(1:A.ptr(end)-1,:);
        A.vals=NC(:,1:A.ptr(end)-1);
        A.tsize=ten.tsize;
        A.tsize(1)=newl;
    else
        l=length(rows);
        temp=zeros(ten.tsize(1),1);
        temp(rows)=1:l;
        [i,j,vals]=tfind(ten);
        i=temp(i);

        if nnz(i)>0
            A=sptensor(i(i~=0),j(i~=0,:),vals(:,i~=0),l,ten.tsize(2:end));
        else
            A=sptensor(l,ten.tsize(2:end),size(ten.vals,1));
        end

        ind=find(i);
    end
end
end

