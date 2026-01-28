function [ coeffs ] = shift_center( coeffs,c0,newc0,order,model )
%The function shifts the center of the approximating power series from c0
%to newc0.
%   coeffs: vector of Polynomial coefficients
%   c0: old center of power series
%   newc0: new center of power series
%   n_f: number of power series
%   n_b: number of terms in power series (excluding symmetric monomials)
%   n_x: number of variables
%   order: Polynomial order of power series
%   U: cell variable with compression matrices
%   W: cell variable with uncompression matrices
%
% © Copyright, Oren Levintal, June 13, 2016.


n_b=model.n_b;
n_x=model.n_x;
U=model.U;
W=model.W;

coeffs=reshape(coeffs,[],n_b);
GH0=coeffs(:,1);
if order==1
    GH1=coeffs(:,2:1+n_x);
    [ GH0,GH1 ] = shiftpoly( newc0,c0,[],GH0,GH1 );
    coeffs=[GH0,GH1];
    coeffs=coeffs(:);
elseif order==2
    GH1=coeffs(:,2:1+n_x);
    GH2=coeffs(:,2+n_x:1+n_x+model.unique2)*W{2};
    [ GH0,GH1,GH2 ] = shiftpoly( newc0,c0,[],GH0,GH1,GH2 );
    coeffs=[GH0,GH1,GH2*U{2}];
    coeffs=coeffs(:);
elseif order==3
    GH1=coeffs(:,2:1+n_x);
    GH2=coeffs(:,2+n_x:1+n_x+model.unique2)*W{2};
    GH3=coeffs(:,2+n_x+model.unique2:1+n_x+model.unique2+model.unique3)*W{3};
    [ GH0,GH1,GH2,GH3 ] = shiftpoly( newc0,c0,[],GH0,GH1,GH2,GH3 );
    coeffs=[GH0,GH1,GH2*U{2},GH3*U{3}];
    coeffs=coeffs(:);
elseif order==4
    GH1=coeffs(:,2:1+n_x);
    GH2=coeffs(:,2+n_x:1+n_x+model.unique2)*W{2};
    GH3=coeffs(:,2+n_x+model.unique2:1+n_x+model.unique2+model.unique3)*W{3};
    GH4=coeffs(:,2+n_x+model.unique2+model.unique3:1+n_x+model.unique2+model.unique3+model.unique4)*W{4};
    [ GH0,GH1,GH2,GH3,GH4 ] = shiftpoly( newc0,c0,[],GH0,GH1,GH2,GH3,GH4 );
    coeffs=[GH0,GH1,GH2*U{2},GH3*U{3},GH4*U{4}];
    coeffs=coeffs(:);
elseif order==5
    GH1=coeffs(:,2:1+n_x);
    GH2=coeffs(:,nchoosek(n_x+1,1)+(1:nchoosek(n_x+1,2)))*W{2};
    GH3=coeffs(:,nchoosek(n_x+2,2)+(1:nchoosek(n_x+2,3)))*W{3};
    GH4=coeffs(:,nchoosek(n_x+3,3)+(1:nchoosek(n_x+3,4)))*W{4};
    GH5=coeffs(:,nchoosek(n_x+4,4)+(1:nchoosek(n_x+4,5)))*W{5};
    [ GH0,GH1,GH2,GH3,GH4,GH5 ] = shiftpoly( newc0,c0,[],GH0,GH1,GH2,GH3,GH4,GH5 );
    coeffs=[GH0,GH1,GH2*U{2},GH3*U{3},GH4*U{4},GH5*U{5}];
    coeffs=coeffs(:);
elseif order==6
    GH1=coeffs(:,2:1+n_x);
    GH2=coeffs(:,nchoosek(n_x+1,1)+(1:nchoosek(n_x+1,2)))*W{2};
    GH3=coeffs(:,nchoosek(n_x+2,2)+(1:nchoosek(n_x+2,3)))*W{3};
    GH4=coeffs(:,nchoosek(n_x+3,3)+(1:nchoosek(n_x+3,4)))*W{4};
    GH5=coeffs(:,nchoosek(n_x+4,4)+(1:nchoosek(n_x+4,5)))*W{5};
    GH6=coeffs(:,nchoosek(n_x+5,5)+(1:nchoosek(n_x+5,6)))*W{6};
    [ GH0,GH1,GH2,GH3,GH4,GH5,GH6 ] = shiftpoly( newc0,c0,[],GH0,GH1,GH2,GH3,GH4,GH5,GH6 );
    coeffs=[GH0,GH1,GH2*U{2},GH3*U{3},GH4*U{4},GH5*U{5},GH6*U{6}];
    coeffs=coeffs(:);
elseif order>=7 
    error('order must not be larger than 6')
end

coeffs=full(coeffs);

end

