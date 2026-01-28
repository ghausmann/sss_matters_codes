function [R0,g,h1,nPhi,nu_vec,model]=eval_model_vec(coeffs,x0,model,params,c0,nep,P,varargin)
% The function evaluates the model residuals of tpvec.
%
% © Copyright, Oren Levintal, October 21, 2016.

order=model.order(1);

ms=0; % no Markov Switching
if ~isempty(varargin)
    coeffsp=varargin{1};
    if ~isempty(coeffsp)
        ms=1; % Markov Switching model
    end
end

if numel(x0)==model.n_x && numel(coeffs)==model.n_theta && numel(params)==model.n_params
    n_s=1;
else % x0,coeffs,params have a row dimension of 1 or n_s
    temp=[size(x0,1),size(coeffs,1),size(params,1),size(c0,1)];
    n_s=max(temp);
    if min(temp)~=n_s
        error('incompatible s dimension')
    end
end

% make sure all inputs are full
coeffsvec=full(reshape(coeffs,n_s,model.n_theta));
x0vec=full(reshape(x0,n_s,model.n_x));

paramsvec=full(reshape(params,n_s,model.n_params));

if model.hybrid==1 || model.hybrid==3
    paramsvec=[zeros(n_s,model.n_e),paramsvec];
end

% eta=full(eta);
if model.hybrid==0 || model.hybrid==1
    etavec=eta_fun(paramsvec,[]);
else
    etavec=zeros(n_s,model.n_x*model.n_e);
end
etavec=reshape(full(etavec),n_s*model.n_x,[]);

c0vec=full(reshape(c0,n_s,model.n_x));
nep=full(nep);
P=full(P);

params=paramsvec;
P=P(:);

n_f=model.n_f; % no. of model conditions
n_x=model.n_x; % no. of state variables
n_y=model.n_y; % no. of control variables
n_x2=model.n_x2; % no. of exogenous state variables
n_x1=model.n_x1; % no. of endogenous state variabels
n_v=model.n_v; % total no. of variables v=(yp,y,xp,x)
n_theta=model.n_theta; % no. of Polynomial coefficients
n_b=model.n_b; % no. of terms in the basis function
n_nodes=size(nep,2); % no. nodes in the discrete shocks
n_z=model.n_z; % no. of auxiliary variables
n_e=model.n_e;

n_prefrows=length(model.prefrows);
n_stochfrows=length(model.stochfrows);

n_params=model.n_params;

if model.hybrid==1 || model.hybrid==3
    n_params=n_params+n_e; % shocks are treated as parameters
end

coeffs=reshape(coeffsvec,n_s,[],n_b);

% current state
nx=x0vec;

nx_c0_vec=x0vec-c0vec;
nx_c0=sptensor(nx_c0_vec(1,:)'); % convert to 1D tensor
nx_c0.vals=nx_c0_vec; % assign the correct values in vals.

if ~isfield(model,'M2')
    model.M2=[];
    model.M3=[];
    model.M3x=[];
    model.TW2=[];
    model.TW3=[];
    model.M4=[];
    model.M4x=[];
    model.M4xx=[];
    model.M4xxx=[];
    if order>=2
        model.M2=fold(spmat2sptensor(model.W{2}*model.W{2}'*model.U{2}'),n_x,n_x);
        model.TW2=fold(spmat2sptensor(model.W{2}),n_x,n_x);
        model.TW2U2=spmat2sptensor(model.W{2}*model.U{2});
    end
    if order>=3
        model.M3=fold(spmat2sptensor(model.W{3}*model.W{3}'*model.U{3}'),n_x,n_x,n_x);
        model.M3x=fold(spmat2sptensor(3*model.W{3}*kron(model.W{2}'*model.U{2}',speye(n_x))),n_x,n_x,n_x);
        model.TW3=fold(spmat2sptensor(model.W{3}),n_x,n_x,n_x);
        model.TW3U3=spmat2sptensor(model.W{3}*model.U{3});
    end
    if order>=4
        model.M4=fold(spmat2sptensor(model.W{4}*model.W{4}'*model.U{4}'),n_x,n_x,n_x,n_x);
        model.M4x=fold(spmat2sptensor(4*model.W{4}*kron(model.W{3}'*model.U{3}',speye(n_x))),n_x,n_x,n_x,n_x);
        model.M4xx=fold(spmat2sptensor(12*model.W{4}*kron(model.W{2}'*model.U{2}',speye(n_x^2))),n_x,n_x,n_x,n_x);
        model.M4xxx=fold(spmat2sptensor(24*model.W{4}),n_x,n_x,n_x,n_x);
    end
end

[X_vecT]=create_X_tensor_no_derivs(order,nx_c0,...
model.M2,model.M3,model.M4,[],model.n_ind,n_s,'vec');

gh_coeffs=reshape(coeffs,n_s,(n_y+n_x1),n_b);
% g_coeffs=reshape(gh_coeffs(:,1:n_y,:),n_s*n_y,n_b);

gh_coeffsT=sptensor(ones(n_y+n_x1,n_b)); 
tempvals=gh_coeffs; %n_s,n_y+n_x1,n_b
tempvals=permute(tempvals,[1,3,2]); % n_s,n_b,n_y+n_x1
tempvals=reshape(tempvals,n_s,n_b*(n_y+n_x1));
gh_coeffsT.vals=tempvals;

[gh1]=contraction1(gh_coeffsT,X_vecT,[],model.n_ind,n_s,'vec');
    
g=takerows(gh1,1:n_y);
h1=takerows(gh1,n_y+(1:n_x1));

% control vars
ny=g.vals;

% expected values: 

% The residual function R is calculated in two steps. The first step
% calculates the nonstochastic (predetermined) rows of f (pref). The second step calculates
% the expected value of the stochastic rows of f (stochf).

% Step 1. predetermined equations (pref)

% preallocate
h=zeros(n_s,n_x);

% build h(x)

h(:,1:n_x1)=h1.vals; % predetermined endogenous state vars

if model.hybrid==0
    nPhi=Phi_fun(nx,params);
    h(:,n_x1+1:end)=nPhi; % expected value of exogenous state vars. shocks are added later
end

% next period state
nxp=h;

% evaluate residuals R0
nv=[zeros(n_s,n_y),ny,nxp,nx]; % all variables with zeros for the stochastic vars.

n_u=model.n_u;
nu=zeros(n_s,n_u);
npreu=preu_fun(nv(:,model.preuvars),params); % all predetermined u
nu(:,model.preurows)=npreu;
nz=[nv,nu];
params(:,model.logical_params_loc)=double(logical(nz(:,model.logical_zvars)>=0));

fname=['pretilf' model.fname '_fun'];
pref=feval(fname,nz(:,model.pretilfzvars),params); % since i already have z, i use pretilf to evaluate f.

preR0=pref;

% Step 2. expected value of stochastic equations (stochf) 

% vectorized expressions
%create stochy: index of future control variables that affect stochastic equations
if isempty(model.stochfrows)
    EstochR0=zeros(n_s,0);
    nu_vec=reshape(repmat(nu(:)',n_nodes,1),n_nodes*n_s,[]);

else
    nx_vec=repmat(nx(:)',n_nodes,1); % n_nodes=the number of nodes of the discrete shocks.
    nx_vec=reshape(nx_vec,n_nodes*n_s,[]);

    % next period state
    % next period state
    if model.hybrid==0
        nxp_vec=repmat(h(:),1,n_nodes)+reshape(etavec,[],n_e)*nep; %n_s*n_x,n_nodes
    elseif model.hybrid==1
        h_vec=repmat(h(:),1,n_nodes); %n_s*n_x,n_nodes
        
        params_vec=[repmat(nep',n_s,1),reshape(repmat(reshape(params(:,n_e+1:end),1,[]),n_nodes,1),n_nodes*n_s,[])]; %n_nodes*n_s,n_params
        uname=['Phi' model.uname '_fun'];
        tempvars=[nx,ny,nu(:,model.preurows)]; %n_s,n_x+n_y
        tempvars=reshape(repmat(tempvars(:)',n_nodes,1),n_nodes*n_s,[]); %n_s*n_nodes,n_x+n_y+n_preurows
        Phi_vec=feval(uname,tempvars,params_vec); %n_s*n_nodes,n_x2
        Phi_vec=permute(reshape(Phi_vec,n_s,n_nodes,n_x2),[1,3,2]); %n_s,n_x2,n_nodes
        nPhi=Phi_vec;
        h_vec=reshape(h_vec,n_s,n_x,n_nodes);
        h_vec(:,n_x1+1:end,:)=Phi_vec; % expected value of exogenous state vars. shocks are added later
        h_vec=reshape(h_vec,n_s*n_x,n_nodes);
        nxp_vec=h_vec+reshape(etavec,[],n_e)*nep; %n_s*n_x,n_nodes
    elseif model.hybrid==3
        h_vec=repmat(h(:),1,n_nodes); %n_s*n_x,n_nodes
        
        params_vec=[repmat(nep',n_s,1),reshape(repmat(reshape(params(:,n_e+1:end),1,[]),n_nodes,1),n_nodes*n_s,[])]; %n_nodes*n_s,n_params
        uname=['Phi' model.uname '_fun'];
        tempvars=nz(:,model.Phizvars); %n_s,n_x+n_y
        tempvars=reshape(repmat(tempvars(:)',n_nodes,1),n_nodes*n_s,[]); %n_nodes*n_s,n_x+n_y+n_preurows
        Phi_vec=feval(uname,tempvars,params_vec); %n_nodes*n_s,n_x2
        Phi_vec=permute(reshape(Phi_vec,n_nodes,n_s,n_x2),[2,3,1]); %n_s,n_x2,n_nodes
        nPhi=Phi_vec;
        h_vec=reshape(h_vec,n_s,n_x,n_nodes);
        h_vec(:,n_x1+1:end,:)=Phi_vec; 
        h_vec=reshape(h_vec,n_s*n_x,n_nodes);
        nxp_vec=h_vec+reshape(etavec,[],n_e)*nep; %n_s*n_x,n_nodes        
    end
    
    nxp_c0_vec=nxp_vec-reshape(repmat(c0vec,1,n_nodes),n_s*n_x,n_nodes);

    nxp_vec_permuted=reshape(permute(reshape(nxp_vec,n_s,n_x,n_nodes),[3,1,2]),n_nodes*n_s,n_x);

    tempxp=1:n_x;
    tempxp=tempxp(:);
    nxp_c0=sptensor(tempxp); % convert to 1D tensor
    tempvals=reshape(nxp_c0_vec,n_s,n_x,n_nodes);
    tempvals=permute(tempvals,[3,1,2]);
    tempvals=reshape(tempvals,n_nodes*n_s,n_x);
    nxp_c0.vals=tempvals; % assign the correct values in vals.


[Xp_vecT]=create_X_tensor_no_derivs(order,nxp_c0,...
    model.M2,model.M3,model.M4,[],model.n_ind,n_s*n_nodes,'vec');

tempv=zeros(n_v,1);
tempv(model.stochfvars)=model.stochfvars;
tempv(n_y+1:end)=0;
stochy=nonzeros(tempv);
n_stochy=length(stochy);
% stochg_coeffs=coeffs(:,stochy,:); 
if ms==0
    stochg_coeffs=coeffs(:,stochy,:); 
else
    coeffspvec=full(reshape(coeffsp,n_s,[],n_b));
    stochg_coeffs=coeffspvec(:,stochy,:); 
end
    stochg_coeffsT=sptensor(ones(n_stochy,n_b)); 
    tempvals=stochg_coeffs; %n_s,n_stochy,n_b
    tempvals=permute(tempvals,[1,3,2]); % n_s,n_b,n_stochy
    tempvals=repmat(tempvals(:)',n_nodes,1); % n_nodes,n_s*n_b*n_stochy
    tempvals=reshape(tempvals,n_nodes*n_s,n_b*n_stochy);
    stochg_coeffsT.vals=tempvals;

    [stochgp_vecT]=contraction1(stochg_coeffsT,Xp_vecT,[],model.n_ind,n_s*n_nodes,'vec');
    
    gp_vec=changerows(stochgp_vecT,stochy,n_y);
    temp=zeros(n_nodes*n_s,n_y);
    temp(:,stochy)=gp_vec.vals;
    gp_vec=temp; % n_nodes*n_s,n_y

% control vars in t+1
nyp_vec=gp_vec; 

% evaluate residuals R0
    nv_vec=[nyp_vec,reshape(repmat(ny(:)',n_nodes,1),n_nodes*n_s,[]),nxp_vec_permuted,nx_vec];

%     nstochv_vec=nv_vec(stochfvars,:); % stochastic vars.

    if model.hybrid==0
        params_vec=reshape(repmat(params(:)',n_nodes,1),n_nodes*n_s,n_params);
    end

    nstochu_vec=stochu_fun(nv_vec(:,model.stochuvars),params_vec); % all stochastic u vars
    nu_vec=reshape(repmat(nu(:)',n_nodes,1),n_nodes*n_s,[]);

    nu_vec(:,model.stochurows)=nstochu_vec;
    nz_vec=[nv_vec,nu_vec];
    params_vec(:,model.logical_params_loc)=double(logical(nz_vec(:,model.logical_zvars)>=0));

    fname=['stochtilf' model.fname '_fun'];
    stochf_vec=feval(fname,nz_vec(:,model.stochtilfzvars),params_vec);


    EstochR0=P'*reshape(stochf_vec,n_nodes,[]); %expected value of stochR0
    EstochR0=reshape(EstochR0,n_s,[]);
end
%     tempT=sptensor(zeros(n_stochfrows,1));
%     tempT.vals=full(EstochR0);
%     EstochR0=tempT;

    R0=zeros(n_s,n_f);
R0(:,model.prefrows)=preR0;
R0(:,model.stochfrows)=EstochR0;

nu_vec=reshape(nu_vec,n_nodes,n_s,[]);
model.coeffs=coeffs;
model.nv=nv;

