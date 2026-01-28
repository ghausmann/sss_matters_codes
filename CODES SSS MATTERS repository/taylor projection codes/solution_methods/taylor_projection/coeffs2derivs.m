function [ derivs0,derivs1,derivs2,derivs3,derivs4,derivs5 ] = coeffs2derivs( coeffs,c0,x0,model )
%The function transforms the vector of unique polynomial coefficients (coeffs) into
%derivatives at x0.

derivs0=[];
derivs1=[];
derivs2=[];
derivs3=[];
derivs4=[];
derivs5=[];

n_b=model.n_b;
n_x=model.n_x;
order=model.order(1);

W=model.W;

coeffs=reshape(coeffs,[],n_b);
newc0=x0;
GH0=coeffs(:,1);
if order==1
    GH1=coeffs(:,2:1+n_x);
    [ GH0,GH1 ] = shiftpoly( newc0,c0,[],GH0,GH1 );
    derivs0=GH0;
    derivs1=GH1;
elseif order==2
    GH1=coeffs(:,2:1+n_x);
    GH2=coeffs(:,2+n_x:1+n_x+model.unique2)*W{2};
    [ GH0,GH1,GH2 ] = shiftpoly( newc0,c0,[],GH0,GH1,GH2 );
    derivs0=GH0;
    derivs1=GH1;
    derivs2=reshape(GH2*2,[],n_x,n_x);
elseif order==3
    GH1=coeffs(:,2:1+n_x);
    GH2=coeffs(:,2+n_x:1+n_x+model.unique2)*W{2};
    GH3=coeffs(:,2+n_x+model.unique2:1+n_x+model.unique2+model.unique3)*W{3};
    [ GH0,GH1,GH2,GH3 ] = shiftpoly( newc0,c0,[],GH0,GH1,GH2,GH3 );
    derivs0=GH0;
    derivs1=GH1;
    derivs2=reshape(GH2*2,[],n_x,n_x);
    derivs3=reshape(GH3*6,[],n_x,n_x,n_x);
elseif order==4
    GH1=coeffs(:,2:1+n_x);
    GH2=coeffs(:,2+n_x:1+n_x+model.unique2)*W{2};
    GH3=coeffs(:,2+n_x+model.unique2:1+n_x+model.unique2+model.unique3)*W{3};
    GH4=coeffs(:,2+n_x+model.unique2+model.unique3:1+n_x+model.unique2+model.unique3+model.unique4)*W{4};
    [ GH0,GH1,GH2,GH3,GH4 ] = shiftpoly( newc0,c0,[],GH0,GH1,GH2,GH3,GH4 );
    derivs0=GH0;
    derivs1=GH1;
    derivs2=reshape(GH2*2,[],n_x,n_x);
    derivs3=reshape(GH3*6,[],n_x,n_x,n_x);
    derivs4=reshape(GH4*24,[],n_x,n_x,n_x,n_x);
elseif order==5
    GH1=coeffs(:,2:1+n_x);
    GH2=coeffs(:,nchoosek(n_x+1,1)+(1:nchoosek(n_x+1,2)))*W{2};
    GH3=coeffs(:,nchoosek(n_x+2,2)+(1:nchoosek(n_x+2,3)))*W{3};
    GH4=coeffs(:,nchoosek(n_x+3,3)+(1:nchoosek(n_x+3,4)))*W{4};
    GH5=coeffs(:,nchoosek(n_x+4,4)+(1:nchoosek(n_x+4,5)))*W{5};
    [ GH0,GH1,GH2,GH3,GH4,GH5 ] = shiftpoly( newc0,c0,[],GH0,GH1,GH2,GH3,GH4,GH5 );
    derivs0=GH0;
    derivs1=GH1;
    derivs2=reshape(GH2*2,[],n_x,n_x);
    derivs3=reshape(GH3*6,[],n_x,n_x,n_x);
    derivs4=reshape(GH4*24,[],n_x,n_x,n_x,n_x);
    derivs5=reshape(GH5*120,[],n_x,n_x,n_x,n_x,n_x);
else 
    error('order must not be larger than 5')
end

end

