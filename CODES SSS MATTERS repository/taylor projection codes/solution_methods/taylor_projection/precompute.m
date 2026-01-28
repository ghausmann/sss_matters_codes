function [ g_theta,h_theta,gx_theta,hx_theta,gxxc_theta,hxxc_theta,gxxxc_theta,hxxxc_theta,gxxxxc_theta,hxxxxc_theta,model ] = precompute(x0,c0,model,order )
%This function precomputese expressions that are independent of coeffs.
%
% © Copyright, Oren Levintal, June 13, 2016.

x0=full(x0);
c0=full(c0);

if ~isequal(x0,c0)
    error('x0 is not equal to c0')
end

n_x=model.n_x;
n_x2=model.n_x2;
n_y=model.n_y;
n_f=model.n_x1+model.n_y; % n_f here is the number of endogenous variables, which might be different than the number of equations in HA models
n_theta=model.n_theta;
U=model.U;
W=model.W;

if ~isfield(model,'unique2')
    model.unique2=[];
end
if ~isfield(model,'unique3')
    model.unique3=[];
end
if ~isfield(model,'unique4')
    model.unique4=[];
end

[X,Xx,Xxx,Xxxx,Xxxxx]=create_X(order,n_x,x0,c0,W{2},model.unique2,W{3},model.unique3,W{4},model.unique4);
X=spmat2sptensor(X); Xxxc=[];Xxxxc=[]; Xxxxxc=[];

if order>=1
    Xx=spmat2sptensor(Xx);
end
if order>=2
    Xxxc=Xxx*U{2};
    Xxxc=spmat2sptensor(Xxxc);
end
if order>=3
    Xxxxc=Xxxx*U{3};
    Xxxxc=spmat2sptensor(Xxxxc);
end
if order>=4
    Xxxxxc=Xxxxx*U{4};
    Xxxxxc=spmat2sptensor(Xxxxxc);
end

if ~isfield(model,'ind')
    model.ind=cell(model.totindi,1);
end
indi=1;

if ~isfield(model,'n_ind')
    n_ind=1;
else
    n_ind=model.n_ind;
end

hx_theta=[]; hxxc_theta=[]; hxxxc_theta=[]; hxxxxc_theta=[];
gx_theta=[]; gxxc_theta=[]; gxxxc_theta=[]; gxxxxc_theta=[];

IfT=spteye(n_f); % this is the identity matrix

[gh_theta,model.ind{indi}]=tkron(ttranspose(X),IfT,model.ind{indi},n_ind);
indi=indi+1;
gh_theta=unfold(gh_theta);
g_theta=takerows(gh_theta,1:n_y);
h_theta=vconcat(takerows(gh_theta,n_y+1:n_f),sptensor(n_x2,n_theta));

if order>=1 
    [ghx_theta,model.ind{indi}]=tkron(ttranspose(Xx),IfT,model.ind{indi},n_ind);
    indi=indi+1;
    ghx_theta=ptr2col(ptr2d(unfold(ghx_theta),n_f,n_x),1); %n_y,n_x,n_theta
    gx_theta=takerows(ghx_theta,1:n_y);
    hx_theta=vconcat(takerows(ghx_theta,n_y+1:n_f),sptensor(n_x2,[n_x,n_theta]));
end    

if order>=2 
    [ghxxc_theta,model.ind{indi}]=tkron(ttranspose(Xxxc),IfT,model.ind{indi},n_ind);
    indi=indi+1;
    ghxxc_theta=ptr2col(ptr2d(unfold(ghxxc_theta),n_f,model.unique2),1); %n_y,unique2,n_theta
    gxxc_theta=takerows(ghxxc_theta,1:n_y);
    hxxc_theta=vconcat(takerows(ghxxc_theta,n_y+1:n_f),sptensor(n_x2,[model.unique2,n_theta]));
end    

if order>=3 
    [ghxxxc_theta,model.ind{indi}]=tkron(ttranspose(Xxxxc),IfT,model.ind{indi},n_ind);
    indi=indi+1;
    ghxxxc_theta=ptr2col(ptr2d(unfold(ghxxxc_theta),n_f,model.unique3),1); %n_y,unique3,n_theta
    gxxxc_theta=takerows(ghxxxc_theta,1:n_y);
    hxxxc_theta=vconcat(takerows(ghxxxc_theta,n_y+1:n_f),sptensor(n_x2,[model.unique3,n_theta]));
end 

if order>=4 
%     [ghxxxxc_theta,model.ind{indi}]=tkron(ttranspose(Xxxxxc),IfT,model.ind{indi},n_ind);
    [ghxxxxc_theta]=tkron(ttranspose(Xxxxxc),IfT);
%     indi=indi+1;
    ghxxxxc_theta=ptr2col(ptr2d(unfold(ghxxxxc_theta),n_f,model.unique4),1); %n_y,unique4,n_theta
    gxxxxc_theta=takerows(ghxxxxc_theta,1:n_y);
    hxxxxc_theta=vconcat(takerows(ghxxxxc_theta,n_y+1:n_f),sptensor(n_x2,[model.unique4,n_theta]));
end 

end

