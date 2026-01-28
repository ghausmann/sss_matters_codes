function [X,Xx,Xxx,Xxxx,Xxxxx,Xxxxxx,ind]=create_X_tensor(N,x_c0,M2,M3,M3x,W2,W3,M4,M4x,M4xx,M4xxx,unique2,unique3,unique4,ind,n_ind,maxload,sum)
%
% © Copyright, Oren Levintal, June 13, 2016.
% changed on January 23, 2017 - extended to 4th order
if isempty(ind)
    ind=cell(9,1);
end
n_x=x_c0.tsize(1);

Xx=[];Xxx=[]; Xxxx=[]; Xxxxx=[]; Xxxxxx=[];
n_s=size(x_c0.vals,1);
temp=sptensor(1);
temp.vals=repmat(temp.vals,n_s,1);
X=temp;
if N==0
    Xx=sptensor(1,n_x,n_s);
elseif N>=1
    X=vconcat(X,x_c0);

    Ix=spteye(n_x,n_s);
    Xx=vconcat(sptensor(1,n_x,n_s),Ix);
end
if N==1
    Xxx=sptensor(1+n_x,n_x^2,n_s);
elseif N>=2
    [temp,ind{1}]=contraction2(M2,x_c0,x_c0,ind{1},n_ind,maxload,sum);
    X=vconcat(X,unfold(temp));
    twoW2=multscalar(W2,2);
    twoW2.vals=repmat(twoW2.vals,n_s,1);
    [tempx,ind{2}]=contraction2(twoW2,x_c0,Ix,ind{2},n_ind,maxload,sum);
    Xx=vconcat(Xx,unfold(tempx));
    Xxx=vconcat(sptensor(1+n_x,n_x^2,n_s),unfold(twoW2));
end
if N==2
    Xxxx=sptensor(1+n_x+unique2,n_x^3,n_s);
elseif N>=3
    [temp,ind{3}]=contraction3(M3,x_c0,x_c0,x_c0,ind{3},n_ind,maxload,sum);
    X=vconcat(X,unfold(temp));
    [tempx,ind{4}]=contraction3(M3x,x_c0,x_c0,Ix,ind{4},n_ind,maxload,sum);
    Xx=vconcat(Xx,unfold(tempx));
    Ix6=multscalar(Ix,6);
    [tempxx,ind{5}]=contraction3(W3,x_c0,Ix,Ix6,ind{5},n_ind,maxload,sum);

    Xxx=vconcat(Xxx,unfold(tempxx));
    sixW3=multscalar(W3,6);
    sixW3.vals=repmat(sixW3.vals,n_s,1);
    Xxxx=vconcat(sptensor(1+n_x+unique2,n_x^3,n_s),unfold(sixW3));
end
if N==3
    Xxxxx=sptensor(1+n_x+unique2+unique3,n_x^4,n_s);
elseif N>=4
    [temp,ind{6}]=contraction4(M4,x_c0,x_c0,x_c0,x_c0,ind{6},n_ind,maxload,sum);
    X=vconcat(X,unfold(temp));
    [tempx,ind{7}]=contraction4(M4x,x_c0,x_c0,x_c0,Ix,ind{7},n_ind,maxload,sum);
    Xx=vconcat(Xx,unfold(tempx));
    [tempxx,ind{8}]=contraction4(M4xx,x_c0,x_c0,Ix,Ix,ind{8},n_ind,maxload,sum);
    Xxx=vconcat(Xxx,unfold(tempxx));
    [tempxxx,ind{9}]=contraction4(M4xxx,x_c0,Ix,Ix,Ix,ind{9},n_ind,maxload,sum);
    Xxxx=vconcat(Xxxx,unfold(tempxxx));
    M4xxx.vals=repmat(M4xxx.vals,n_s,1);
    Xxxxx=vconcat(sptensor(1+n_x+unique2+unique3,n_x^4,n_s),unfold(M4xxx));
end
if N==4
    Xxxxxx=sptensor(1+n_x+unique2+unique3+unique4,n_x^5,n_s);
end
