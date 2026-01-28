function model=prepare_tpvec(f,Phi,yp,y,xp,x,symparams,eta,order,varargin)

% This file prepares codes and data that are used later to compute the
% nonlinear system of Taylor projection and the Jacobian.
% Input args:
%   f = symbolic vector of the model conditions (Ef=0)
%   Phi = symbolic vector of the lower block of h(x). This is the evolution law
%   of the exogenoust state variables, or any other state variables with a
%   known policy function.
%   yp,y,xp,x = symbolic vectors of the control and state varaibles for
%   current (y,x) and future (yp,xp) periods.
%   symparams = symbolic vector of the model parameters.
%   eta = the matrix eta as in Schmitt-Grohe and Uribe (2004).
%   order = order of the approximiating Polynomials. Can be 1, 2 or 3.
% Optional arguments:
%   subsfuns,subsvars = subsvars (the second optional argument) is a symbolic
%   vector of auxiliary variables that are eventually substituted out.
%   subsfuns (the first optional argument) is a symbolic vector of the
%   expressions of subsvars (see example in the documentation). Use
%   subsfuns and subsvars to speed up differentiation.
% Output arg:
%   model = a structure variable with data, that is used to compute the
%   system.
%
% © Copyright, Oren Levintal, June 13, 2016.
% changed on March 19, 2017 - the function is vectorized and extended to
% fourth order.



model.n_e=size(eta,2);

f=sym(f(:));
Phi=sym(Phi(:));
eta=sym(eta);

x=sym(x(:));
xp=sym(xp(:));
y=sym(y(:));
yp=sym(yp(:));

if length(x)~=length(xp)
    error('wrong size of current or future state variables')
end
if length(y)~=length(yp)
    error('wrong size of current or future control variables')
end

    
% SIZE VARIABLES
n_f=length(f); n_x=length(x); n_y=length(y); n_x2=length(Phi); n_x1=n_x-n_x2;
v=[yp(:); y(:); xp(:); x(:)];
n_v=length(v);

if size(eta,1)~=n_x
    error('wrong number of rows for matrix eta')
end
if nnz(eta(1:n_x1,:))~=0
    error('the first n_x1 rows of eta should be zero')
end

% if n_f>n_y+n_x1
%     error('number of model conditions exceeds the number of endogenous variables');
% elseif n_f<n_y+n_x1
%     error('number of endogenous variables exceeds the number of model conditions');
% end

% if order is a vector of length 2, the first element is the order of
% Taylor projection and the second is the order of perturbation.
if length(order)==2
    pert_order=order(2);
    if pert_order>4
        error('perturbation order cannot exceed 4');
    end
    order=order(1);
else
    pert_order=order;
end
model.order=[order,pert_order];

n_b=1;
for k=1:order
    n_b=n_b+nchoosek(n_x+k-1,k); % number of parameters in the basis function
end

n_theta=n_b*(n_y+n_x1); % total number of parameters to be solved

% if order>=2
    model.unique2=nchoosek(n_x+1,2);
% end
% if order>=3
    model.unique3=nchoosek(n_x+2,3);
% end
% if order>=4
    model.unique4=nchoosek(n_x+3,4);
% end

% fname

if nargin>11
    fname=varargin{3};
    if ~ischar(fname)
        error('function name must be a string')
    end
else
    fname='';
end
if nargin>12
    uname=varargin{4};
    if ~ischar(fname)
        error('function name must be a string')
    end
else
    uname='';
end
model.fname=fname;
model.uname=uname;

if nargin>13
    logicalparams=varargin{5};
    logicalvars=varargin{6};
    if length(logicalparams)~=length(logicalvars)
        error('wrong size of logical parameters or variables')
    end
else
    logicalparams=sym([]);
    logicalvars=sym([]);
end

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
    if length(pi_)~=length(u)
        error('wrong size of substituted variables or functions')
    end
    if ~isequal(sort(u),sort(unique(u)))
        error('substituted variables are not uniquely determined')
    end
end

% check that auxiliary functions do not depend on logical parameters
if ~isempty(logicalparams) && ~isempty(pi_)
    if ~isequal(pi_,subs(pi_,logicalparams,sym(zeros(size(logicalparams)))))
        error('auxiliary variables must not depend on logical parameters')
    end
end

% check that Phi function does not depend on logical parameters
if ~isempty(logicalparams) && ~isempty(Phi)
    if ~isequal(Phi,subs(Phi,logicalparams,sym(zeros(size(logicalparams)))))
        error('function Phi must not depend on logical parameters')
    end
end

% check that all symbolic variables/parameters are defined

actualsymvars=symvar([f;pi_;eta(:);Phi(:)]);
definedvars=[v(:);u(:);symparams(:);logicalparams(:)];
actualsymvars=sort(actualsymvars(:));
definedvars=sort(definedvars);

n_definedvars=length(definedvars);
j=0;
for i=1:length(actualsymvars)
    j=j+1;
    if j>n_definedvars
        error([char(actualsymvars(i)) ' is not defined'])
    end
    while ~isequal(actualsymvars(i),definedvars(j))
        j=j+1;
        if j>n_definedvars
            error([char(actualsymvars(i)) ' is not defined'])
        end
    end
end

% check that logical variables are defined
actualsymvars=unique(logicalvars);
definedvars=[v(:);u(:)];
actualsymvars=sort(actualsymvars(:));
definedvars=sort(definedvars);
n_definedvars=length(definedvars);

j=0;
for i=1:length(actualsymvars)
    j=j+1;
    if j>n_definedvars
        error([char(actualsymvars(i)) ' is not defined'])
    end
    while ~isequal(actualsymvars(i),definedvars(j))
        j=j+1;
        if j>n_definedvars
            error([char(actualsymvars(i)) ' is not defined'])
        end
    end
end


% for i=1:length(logicalvars)
%     if isempty(v(v==logicalvars(i))) && isempty(u(u==logicalvars(i))) 
%         error([char(logicalvars(i)) ' is not defined'])
%     end
% end

symparams=[symparams(:);logicalparams(:)];
model.n_params=numel(symparams);
model.logical_params_loc=model.n_params-length(logicalparams)+1:model.n_params;
model.n_logicalparams=length(logicalparams);

n_u=length(u);
model.n_u=n_u;
tilf=f; % tilf(v,u) is a function of v and u
% n=find_n(pi_,u); % number of substitutions needed to eliminate u
[n,piu0]=find_n_new(pi_,u);
model.uu0=piu0;

% Identify variables that affect u
u_v=u; %u_v which is u as a function of v only u(v)
uu=eye(n_u); % matrix to store which of u affect u
uv=zeros(n_u,n_v);

% piu0=double(logical(jacobian(pi_,u)~=0));
piv0=double(logical(jacobian(pi_,v)~=0));
piu=piu0;
piv=piv0;

for k=1:n 
    uu=uu+piu;
    uv=uv+piv;
    piv=piu*piv0;
    piu=piu*piu0;
    u_v=subs(u_v,u,pi_); % substitute pi_ into itself n times
end

% for k=1:n 
%     u_v=subs(u_v,u,pi_); % substitute pi_ into itself n times
%     uu=uu+logical(jacobian(u_v,u)~=0); % uu is a n_u-by-n_u matrix. the ij element is zero only if ui is independent of uj through all substitutions
%     uv=uv+logical(jacobian(u_v,v)~=0); % uv is a n_u-by-n_v matrix. the ij element is zero only if ui is independent of vj through all substitutions
% end

% Identify stochasic and predetermined functions and variables
if n_x2>0
    stochexog=find(sum(1-logical(eta(n_x1+1:end,:)==0),2)); % stochastic exogenous variables
else
    stochexog=[];
end
preexog=1:n_x2;
preexog(stochexog)=[]; % predetermined exogenous variables

stochvars=[1:n_y,2*n_y+n_x1+stochexog']; % all stochastic variables

tilfv=double(logical(jacobian(tilf,v)~=0));
tilfu=double(logical(jacobian(tilf,u)~=0));
fv=double(logical(tilfv+tilfu*uv~=0)); % logical Jacobian of f w.r.t v 

stochfrows=find(sum( 1-logical(fv(:,stochvars)==0),2)); % stochastic rows of f
prefrows=1:n_f;
prefrows(stochfrows)=[]; % predetermined rows of f

% variables that affect the predetermined and stochastic rows of f

prefvars=find(sum(1-logical(fv(prefrows,:)==0),1));
stochfvars=find(sum(1-logical(fv(stochfrows,:)==0),1));


currentFolder=pwd;
if ~exist([pwd '\files'], 'dir')
   mkdir('files');
end

cd 'files'

% generate a function that computes eta
gen_fun_vec(eta(:),[],symparams,['eta' uname],'row');


% Create an m file for Phi and its derivatives w.r.t x
gen_fun_vec(Phi,symparams,x,['Phi' uname],'row');
if order>=1
    model.Phi_indc=getderivs_c(Phi,x,max(order,pert_order),symparams,'Phi');
end

% Create an m file for u
% newuv=jacobian(u_v,v); % calculate uv again
newuv=uv;
stochurows=find(sum( 1-logical(newuv(:,stochvars)==0),2)); % stochastic rows of u
model.stochurows=stochurows;
preurows=1:n_u;
preurows(stochurows)=[]; % predetermined rows of u
model.preurows=preurows;
preuvars=find(sum( logical(newuv(preurows,:)~=0),1)); % variables of predetermined rows of u
stochuvars=find(sum( logical(newuv(stochurows,:)~=0),1)); % variables of stochastic rows of u

% gen_fun_vec(u_v(preurows),symparams,v(preuvars),['preu' uname],'row');
% gen_fun_vec(u_v(stochurows),symparams,v(stochuvars),['stochu' uname],'row');
preurows_allsubs=logical(sum(uu(preurows,:),1));
gen_fun_vec_subs(u(preurows),symparams,v(preuvars),pi_(preurows_allsubs),u(preurows_allsubs),['preu' uname],'row',piu0(preurows_allsubs,preurows_allsubs));

stochurows_allsubs=logical(sum(uu(stochurows,:),1));
gen_fun_vec_subs(u(stochurows),symparams,v(stochuvars),pi_(stochurows_allsubs),u(stochurows_allsubs),['stochu' uname],'row',piu0(stochurows_allsubs,stochurows_allsubs));

model.preuvars=preuvars;
model.stochuvars=stochuvars;

% build z and find the Jacobian of z
z=[v;u];

% locate the logical variables in z

logical_zvars=zeros(length(logicalvars),1);
for i=1:length(logical_zvars)
   logical_zvars(i)=find(logical(z==logicalvars(i)));
end
model.logical_zvars=logical_zvars;

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

% tilfz=jacobian(tilf,z);
tilfz=[tilfv,tilfu];

% build tilfz
% tilfz=double(logical(tilfz~=0)); % only direct effects
for i=1:n_f
    temp=tilfz(i,:);
    temp(temp~=0)=1:nnz(temp);
    tilfz(i,:)=temp;
end
model.tilfz=tilfz;
model.maxtilfz=intarray(max(tilfz(:)));

% Create an m file for tilf
pretilfzvars=find(sum(1-logical(tilfz(prefrows,:)==0),1));
model.pretilfzvars=pretilfzvars; % variables in z that affect pretilf
gen_fun_vec(tilf(prefrows),symparams,z(pretilfzvars),['pretilf' fname],'row');

stochtilfzvars=find(sum(1-logical(tilfz(stochfrows,:)==0),1));
model.stochtilfzvars=stochtilfzvars; % variables in z that affect stochtilf
gen_fun_vec(tilf(stochfrows),symparams,z(stochtilfzvars),['stochtilf' fname],'row');

% Differentiate PI and create m files

PI=[v;pi_];
prePI=PI;
prePI(n_v+stochurows)=0; % stochastic rows of PI are set to zero. Note: the first n_v rows include stochastic variables but their derivatives w.r.t are ones or zeros.
stochPI=PI;
stochPI([1:n_v,n_v+preurows])=0; % nonstochastic rows are set to zero. See the note above.

[model.prePI_ind_u]=getderivs_tensor(prePI,z,max(order+1,pert_order),symparams,['prePI' uname]); 
[model.stochPI_ind_u]=getderivs_tensor(stochPI,z,max(order+1,pert_order),symparams,['stochPI' uname]);

prefuvars=find(sum(uu(logical(sum(logical(tilfu(prefrows,:)~=0),1)~=0),:),1));
stochfuvars=find(sum(uu(logical(sum(logical(tilfu(stochfrows,:)~=0),1)~=0),:),1));
model.prefuvars=prefuvars;
model.stochfuvars=stochfuvars;

prefzvars=[prefvars,n_v+prefuvars];
stochfzvars=[stochfvars,n_v+stochfuvars];
model.n_prefzvars=length(prefzvars);
model.n_stochfzvars=length(stochfzvars);

[model.stochtilf_ind_u,model.stoch_n]=gen_chainderivs_tensor(tilf(stochfrows),v(stochfvars),u(stochfuvars),pi_(stochfuvars),max(order+1,pert_order),symparams,['stochtilf' fname]);
[model.pretilf_ind_u,model.pre_n]=gen_chainderivs_tensor(tilf(prefrows),v(prefvars),u(prefuvars),pi_(prefuvars),max(order+1,pert_order),symparams,['pretilf' fname]);


% Store variables in struct model
model.prefrows=prefrows;
model.stochfrows=stochfrows;
model.prefvars=prefvars;
model.stochfvars=stochfvars;
model.stochexog=stochexog;
model.preexog=preexog;
model.n_stochexog=length(stochexog);
model.n_preexog=length(preexog);
model.stochfzvars=stochfzvars;
model.prefzvars=prefzvars;

model.n_theta=n_theta;
model.n_b=n_b;

model.n_f=n_f; 
model.n_x=n_x; 
model.n_y=n_y; 
model.n_x1=n_x1;
model.n_x2=n_x2; 

model.n_v=n_v;
model.n_z=n_z;
model.n_u=n_u;

model.n_ind=1;

model.U=cell(4,1);
model.W=cell(4,1);
for i=2:order+1
    [model.U{i},model.W{i}]=create_UW(n_x,i);
end

model.stochzz=zz(model.stochfzvars,model.stochfzvars);
model.prezz=zz(model.prefzvars,model.prefzvars);

model.stochtilfz=tilfz(model.stochfrows,model.stochfzvars);
model.pretilfz=tilfz(model.prefrows,model.prefzvars);

% calculate some big matrices now

if order>=3
    save('tempfile')
    n_prefvars=length(prefvars);
    clearvars -except n_prefvars

    if n_prefvars>0
        [ tempM ] = chainsM( n_prefvars,4 );
    else
        tempM=cell(4,1);
    end
    load('tempfile')
    model.prefvars_chain4c_M2=tempM{2};
    model.prefvars_chain4c_M3=tempM{3};
    model.prefvars_chain4c_M4=tempM{4};
    clear tempM
    save('tempfile')
    
    
    
    n_stochfvars=length(stochfvars);
    clearvars -except n_stochfvars

    if n_stochfvars>0
        [ tempM ] = chainsM( n_stochfvars,4 );
    else
        tempM=cell(4,1);
    end
    load('tempfile')
    model.stochfvars_chain4c_M2=tempM{2};
    model.stochfvars_chain4c_M3=tempM{3};
    model.stochfvars_chain4c_M4=tempM{4};
    clear tempM
    delete('tempfile.mat')
end

if order>=4
    save('tempfile')
    n_prefvars=length(prefvars);
    clearvars -except n_prefvars

    if n_prefvars>0
        [ tempM ] = chainsM( n_prefvars,5 );
    else
        tempM=cell(6,1);
    end
    load('tempfile')
    model.prefvars_chain5c_M1=tempM{1};
    model.prefvars_chain5c_M2=tempM{2};
    model.prefvars_chain5c_M3=tempM{3};
    model.prefvars_chain5c_M4=tempM{4};
    model.prefvars_chain5c_M5=tempM{5};
    model.prefvars_chain5c_M6=tempM{6};

    clear tempM
    save('tempfile')
    
    
    
    n_stochfvars=length(stochfvars);
    clearvars -except n_stochfvars

    if n_stochfvars>0
        [ tempM ] = chainsM( n_stochfvars,5 );
    else
        tempM=cell(6,1);
    end
    load('tempfile')
    model.stochfvars_chain5c_M1=tempM{1};
    model.stochfvars_chain5c_M2=tempM{2};
    model.stochfvars_chain5c_M3=tempM{3};
    model.stochfvars_chain5c_M4=tempM{4};
    model.stochfvars_chain5c_M5=tempM{5};
    model.stochfvars_chain5c_M6=tempM{6};

    clear tempM
    delete('tempfile.mat')
end

totindi=5+8;
totindi=totindi+(model.pre_n-1)*(1+3+4+6)+1+4+6+9;
totindi=totindi+3;
totindi=totindi+4;
totindi=totindi+2;
totindi=totindi+2;
totindi=totindi+1;
totindi=totindi+9;
totindi=totindi+8;
totindi=totindi+8+(model.stoch_n-1)*(1+3+4+6)+1+4+6+9;
totindi=totindi+2;
totindi=totindi+9;
totindi=totindi+7;

model.totindi=totindi;

model.count_pre_n=0;
model.count_stoch_n=0;

rehash;


cd(currentFolder)