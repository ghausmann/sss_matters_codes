function [X,Xx,Xxx,Xxxx,Xxxxx,Xxxxxx]=create_X(N,n_x,x,c0,W2,unique2,W3,unique3,W4,unique4)
% compute the basis function and its derivatives
%
% © Copyright, Oren Levintal, June 13, 2016.
% extended to 4th order on January 19, 2017.

Xx=[];Xxx=[]; Xxxx=[]; Xxxxx=[]; Xxxxxx=[];
if N==0
    X=1;
    Xx=sparse(1,n_x);
elseif N>=1
    x_c0=sparse(x-c0);
    X=[1;x_c0];
    Xx=[sparse(1,n_x);speye(n_x)];
end
if N==1
    Xxx=sparse(1+n_x,n_x^2);
elseif N>=2
    tempx=reshape(reshape(W2,unique2*n_x,n_x)*x_c0,unique2,n_x);
    temp=tempx*x_c0;
    X=[X;temp];
    Xx=[Xx;2*tempx];
    Xxx=[sparse(1+n_x,n_x^2);2*W2];
end
if N==2
    Xxxx=sparse(1+n_x+unique2,n_x^3);
elseif N>=3
    tempxx=reshape(reshape(W3,unique3*n_x^2,n_x)*x_c0,unique3,n_x^2);
    tempx=reshape(reshape(tempxx,unique3*n_x,n_x)*x_c0,unique3,n_x);
    temp=tempx*x_c0;
    X=[X;temp];
    Xx=[Xx;3*tempx];
    Xxx=[Xxx;6*tempxx];
    Xxxx=[sparse(1+n_x+unique2,n_x^3);6*W3];
end
if N==3
    Xxxxx=sparse(1+n_x+unique2+unique3,n_x^4);
elseif N>=4
    tempxxx=reshape(reshape(W4,unique4*n_x^3,n_x)*x_c0,unique4,n_x^3);
    tempxx=reshape(reshape(tempxxx,unique4*n_x^2,n_x)*x_c0,unique4,n_x^2);
    tempx=reshape(reshape(tempxx,unique4*n_x,n_x)*x_c0,unique4,n_x);
    temp=tempx*x_c0;
    X=[X;temp];
    Xx=[Xx;4*tempx];
    Xxx=[Xxx;12*tempxx];
    Xxxx=[Xxxx;24*tempxxx];
    Xxxxx=[sparse(1+n_x+unique2+unique3,n_x^4);24*W4];
end
if N==4
    Xxxxxx=sparse(1+n_x+unique2+unique3+unique4,n_x^5);
end
end
