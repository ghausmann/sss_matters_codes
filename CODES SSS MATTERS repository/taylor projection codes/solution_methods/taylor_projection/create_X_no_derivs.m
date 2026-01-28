function [X]=create_X_no_derivs(N,n_x,x,c0,W2,unique2,W3,unique3,W4,unique4)
% compute the basis function and its derivatives
%
% © Copyright, Oren Levintal, April 25, 2017.


x_c0=sparse(x-c0);
X=[1;x_c0];
if N>=2
    tempx=reshape(reshape(W2,unique2*n_x,n_x)*x_c0,unique2,n_x);
    temp=tempx*x_c0;
    X=[X;temp];
end
if N>=3
    tempxx=reshape(reshape(W3,unique3*n_x^2,n_x)*x_c0,unique3,n_x^2);
    tempx=reshape(reshape(tempxx,unique3*n_x,n_x)*x_c0,unique3,n_x);
    temp=tempx*x_c0;
    X=[X;temp];
end
if N>=4
    tempxxx=reshape(reshape(W4,unique4*n_x^3,n_x)*x_c0,unique4,n_x^3);
    tempxx=reshape(reshape(tempxxx,unique4*n_x^2,n_x)*x_c0,unique4,n_x^2);
    tempx=reshape(reshape(tempxx,unique4*n_x,n_x)*x_c0,unique4,n_x);
    temp=tempx*x_c0;
    X=[X;temp];
end

end
