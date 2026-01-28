function model=diff_fun(f,x,symparams,order,varargin)

% This file differentiates the symbolic function f w.r.t
% the symbolic vector x up to order and returns data used to compute those
% derivatives.
%
% © Copyright, Oren Levintal, September 8, 2017.

model.n_params=numel(symparams);

f=sym(f(:));

x=sym(x(:));

  
% SIZE VARIABLES
n_f=length(f); n_x=length(x);
v=x(:);
n_v=length(v);

if order>5
    error('differentiation is allowed up to fifth order')
end
model.order=order;

% fname

if nargin>4
    fname=varargin{3};
    if ~ischar(fname)
        error('function name must be a string')
    end
else
    fname='';
end

model.fname=fname;

% DEFINE pi_ and u
if isempty(varargin)
    pi_=sym([]);
    u=sym([]);
elseif isempty(varargin{1})
    warning('no substitutions assumed')
    pi_=sym([]);
    u=sym([]);
else
    pi_=varargin{1}; % pi_ is the function pi(v,u)
    pi_=pi_(:);
    u=varargin{2}; % auxiliary variables that are eventually substituted out
    u=u(:);
    if ~isequal(sort(u),unique(u))
        error('substituted variables are not uniquely determined')
    end
end

% check that all symbolic variables/parameters are defined
actualsymvars=symvar([f;pi_]);
for i=1:length(actualsymvars)
    if isempty(v(v==actualsymvars(i))) && isempty(u(u==actualsymvars(i))) && isempty(symparams(symparams==actualsymvars(i)))
        error([char(actualsymvars(i)) ' is not defined'])
    end
end

n_u=length(u);
model.n_u=n_u;
tilf=f; % tilf(v,u) is a function of v and u
n=find_n(pi_,u); % number of substitutions needed to eliminate u

% Identify variables that affect u
u_v=u; %u_v which is u as a function of v only u(v)
uu=eye(n_u); % matrix to store which of u affect u
uv=zeros(n_u,n_v);
for k=1:n 
    u_v=subs(u_v,u,pi_); % substitute pi_ into itself n times
    uu=uu+logical(jacobian(u_v,u)~=0); % uu is a n_u-by-n_u matrix. the ij element is zero only if ui is independent of uj through all substitutions
    uv=uv+logical(jacobian(u_v,v)~=0); % uv is a n_u-by-n_v matrix. the ij element is zero only if ui is independent of vj through all substitutions
end

fv=logical(logical(jacobian(tilf,v)~=0)+logical(jacobian(tilf,u)~=0)*uv~=0); % logical Jacobian of f w.r.t v 
fvars=find(sum(1-logical(fv==0),1));

currentFolder=pwd;
if ~exist([pwd '\fun_' fname], 'dir')
   mkdir(['fun_' fname]);
end

cd(['fun_' fname])

% Create an m file for u
gen_fun_vec(u_v,symparams,v,'u','row');

% build z and find the Jacobian of z
z=[v;u];
n_z=n_v+n_u;
zz=[eye(n_v),zeros(n_v,n_u);uv,uu];
for i=1:n_z
    temp=zz(i,:);
    temp(temp~=0)=1:nnz(temp);
    zz(i,:)=temp;
end
zz=intarray(zz);
model.maxzz=intarray(max(zz(:)));

% build fv similar to zz
fv=double(fv);
for i=1:n_f
    temp=fv(i,:);
    temp(temp~=0)=1:nnz(temp);
    fv(i,:)=temp;
end
model.fv=fv;

tilfz=jacobian(tilf,z);

% build tilfz
tilfz=double(logical(tilfz~=0)); % only direct effects
for i=1:n_f
    temp=tilfz(i,:);
    temp(temp~=0)=1:nnz(temp);
    tilfz(i,:)=temp;
end
model.tilfz=tilfz;
model.maxtilfz=intarray(max(tilfz(:)));

% Create an m file for tilf
tilfzvars=find(sum(1-logical(tilfz==0),1));
model.tilfzvars=tilfzvars; % variables in z that affect tilf
gen_fun_vec(tilf,symparams,z(tilfzvars),['tilf' fname],'row');

% Differentiate PI and create m files

PI=[v;pi_];

[model.PI_ind_u]=getderivs_tensor(PI,z,order,symparams,'PI'); 

fuvars=find(sum(uu(logical(sum(logical(jacobian(tilf,u)~=0),1)~=0),:),1));
model.fuvars=fuvars;

fzvars=[fvars,n_v+fuvars];
model.n_fzvars=length(fzvars);

[model.tilf_ind_u,model.pre_n]=gen_chainderivs_tensor(tilf,v,u,pi_,order,symparams,['tilf' fname]);

% Store variables in struct model
model.fvars=fvars;
model.fzvars=fzvars;

model.n_f=n_f; 
model.n_x=n_x; 

model.n_v=n_v;
model.n_z=n_z;
model.n_u=n_u;

model.n_ind=1;

% calculate some big matrices now

if order>=3
    save('tempfile')
    n_fvars=length(fvars);
    clearvars -except n_fvars

    if n_fvars>0
        [ tempM ] = chainsM( n_fvars,4 );
    else
        tempM=cell(4,1);
    end
    load('tempfile')
    model.fvars_chain4c_M2=tempM{2};
    model.fvars_chain4c_M3=tempM{3};
    model.fvars_chain4c_M4=tempM{4};
    clear tempM
    delete('tempfile.mat')
end

if order>=4
    save('tempfile')
    n_fvars=length(fvars);
    clearvars -except n_fvars

    if n_fvars>0
        [ tempM ] = chainsM( n_efvars,5 );
    else
        tempM=cell(6,1);
    end
    load('tempfile')
    model.fvars_chain5c_M1=tempM{1};
    model.fvars_chain5c_M2=tempM{2};
    model.fvars_chain5c_M3=tempM{3};
    model.fvars_chain5c_M4=tempM{4};
    model.fvars_chain5c_M5=tempM{5};
    model.fvars_chain5c_M6=tempM{6};

    clear tempM

    delete('tempfile.mat')
end

rehash;


cd(currentFolder)