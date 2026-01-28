function model=diff_fun2(f,x,symparams,order,varargin)

model.n_params=numel(symparams);

f=sym(f(:));
x=sym(x(:));

% SIZE VARIABLES
n_f=length(f); n_x=length(x); 
v=[x(:)];
n_v=length(v);

model.order=[order];

% fname

if nargin==7
    fname=varargin{3};
    if ~ischar(fname)
        error('function name must be a string')
    end
else
    fname='new';
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
    if length(pi_)~=length(u)
        error('wrong size of substituted variables or functions')
    end
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
% [n,u_v]=find_n(pi_,u); % number of substitutions needed to eliminate u
n=find_n_new(pi_,u);

% Identify variables that affect u
uu=eye(n_u); % matrix to store which of u affect u
uv=zeros(n_u,n_v);

piu0=double(logical(jacobian(pi_,u)~=0));
piv0=double(logical(jacobian(pi_,v)~=0));
piu=piu0;
piv=piv0;

for k=1:n 
    uu=uu+piu;
    uv=uv+piv;
    piv=piu*piv0;
    piu=piu*piu0;
end
% uu_=uu;
% uv_=uv;
% % older version
% u_v=u; %u_v which is u as a function of v only u(v)
% uu=eye(n_u); % matrix to store which of u affect u
% uv=zeros(n_u,n_v);
% for k=1:n 
%     u_v=subs(u_v,u,pi_); % substitute pi_ into itself n times
%     uu=uu+logical(jacobian(u_v,u)~=0); % uu is a n_u-by-n_u matrix. the ij element is zero only if ui is independent of uj through all substitutions
%     uv=uv+logical(jacobian(u_v,v)~=0); % uv is a n_u-by-n_v matrix. the ij element is zero only if ui is independent of vj through all substitutions
% end

% Identify stochasic and predetermined functions and variables
% if n_x2>0
%     stochexog=find(sum(1-logical(eta(n_x1+1:end,:)==0),2)); % stochastic exogenous variables
% else
%     stochexog=[];
% end
% preexog=1:n_x2;
% preexog(stochexog)=[]; % predetermined exogenous variables
% 
% stochvars=[1:n_y,2*n_y+n_x1+stochexog']; % all stochastic variables
fv=ones(n_f,n_v);
% fv=logical(logical(jacobian(tilf,v)~=0)+logical(jacobian(tilf,u)~=0)*uv~=0); % logical Jacobian of f w.r.t v 

% stochfrows=find(sum( 1-logical(fv(:,stochvars)==0),2)); % stochastic rows of f
prefrows=1:n_f;
% prefrows(stochfrows)=[]; % predetermined rows of f

% variables that affect the predetermined and stochastic rows of f

prefvars=find(sum(1-logical(fv(prefrows,:)==0),1));
% stochfvars=find(sum(1-logical(fv(stochfrows,:)==0),1));


currentFolder=pwd;
if ~exist([pwd '\files'], 'dir')
   mkdir('files');
end

cd 'files'

% generate a function that computes eta
% gen_fun_vec(eta(:),[],symparams,'eta','row');


% Create an m file for Phi and its derivatives w.r.t x
% gen_fun_vec(Phi,symparams,x,'Phi','row');
% if order>=1
%     model.Phi_indc=getderivs_c(Phi,x,max(order,pert_order),symparams,'Phi');
% end

% Create an m file for u
% newuv=jacobian(u_v,v); % calculate uv again
newuv=uv;
% stochurows=find(sum( 1-logical(newuv(:,stochvars)==0),2)); % stochastic rows of u
% model.stochurows=stochurows;
preurows=1:n_u;
% preurows(stochurows)=[]; % predetermined rows of u
model.preurows=preurows;
preuvars=find(sum( logical(newuv(preurows,:)~=0),1)); % variables of predetermined rows of u
% stochuvars=find(sum( logical(newuv(stochurows,:)~=0),1)); % variables of stochastic rows of u
gen_fun_vec_subs(u(preurows),symparams,v(preuvars),pi_,u,[fname 'u'],'row');
% gen_fun_vec(u_v(stochurows),symparams,v(stochuvars),'stochu','row');
model.preuvars=preuvars;
% model.stochuvars=stochuvars;

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
pretilfzvars=find(sum(1-logical(tilfz(prefrows,:)==0),1));
model.pretilfzvars=pretilfzvars; % variables in z that affect pretilf
gen_fun_vec(tilf(prefrows),symparams,z(pretilfzvars),[fname 'f'],'row');

% stochtilfzvars=find(sum(1-logical(tilfz(stochfrows,:)==0),1));
% model.stochtilfzvars=stochtilfzvars; % variables in z that affect stochtilf
% gen_fun_vec(tilf(stochfrows),symparams,z(stochtilfzvars),['stochtilf' fname],'row');

% Differentiate PI and create m files

PI=[v;pi_];
prePI=PI;
% prePI(n_v+stochurows)=0; % stochastic rows of PI are set to zero. Note: the first n_v rows include stochastic variables but their derivatives w.r.t are ones or zeros.
% stochPI=PI;
% stochPI([1:n_v,n_v+preurows])=0; % nonstochastic rows are set to zero. See the note above.

[model.prePI_ind_u]=getderivs_tensor(prePI,z,order,symparams,[fname 'PI']); 
% [model.stochPI_ind_u]=getderivs_tensor(stochPI,z,max(order+1,pert_order),symparams,'stochPI');

prefuvars=find(sum(uu(logical(sum(logical(jacobian(tilf(prefrows),u)~=0),1)~=0),:),1));
% stochfuvars=find(sum(uu(logical(sum(logical(jacobian(tilf(stochfrows),u)~=0),1)~=0),:),1));
model.prefuvars=prefuvars;
% model.stochfuvars=stochfuvars;

prefzvars=[prefvars,n_v+prefuvars];
% stochfzvars=[stochfvars,n_v+stochfuvars];
model.n_prefzvars=length(prefzvars);
% model.n_stochfzvars=length(stochfzvars);

% [model.stochtilf_ind_u,model.stoch_n]=gen_chainderivs_tensor(tilf(stochfrows),v(stochfvars),u(stochfuvars),pi_(stochfuvars),max(order+1,pert_order),symparams,['stochtilf' fname]);
[model.pretilf_ind_u,model.pre_n]=gen_chainderivs_tensor(tilf(prefrows),v(prefvars),u(prefuvars),pi_(prefuvars),order,symparams,[fname 'tilf']);


% Store variables in struct model
model.prefrows=prefrows;
% model.stochfrows=stochfrows;
model.prefvars=prefvars;
% model.stochfvars=stochfvars;
% model.stochexog=stochexog;
% model.preexog=preexog;
% model.n_stochexog=length(stochexog);
% model.n_preexog=length(preexog);
% model.stochfzvars=stochfzvars;
model.prefzvars=prefzvars;

% model.n_theta=n_theta;
% model.n_b=n_b;

model.n_f=n_f; 
model.n_x=n_x; 
% model.n_y=n_y; 
% model.n_x1=n_x1;
% model.n_x2=n_x2; 

model.n_v=n_v;
model.n_z=n_z;
model.n_u=n_u;

model.n_ind=1;

% model.U=cell(4,1);
% model.W=cell(4,1);
% for i=2:order+1
%     [model.U{i},model.W{i}]=create_UW(n_x,i);
% end

% model.stochzz=zz(model.stochfzvars,model.stochfzvars);
model.prezz=zz(model.prefzvars,model.prefzvars);

% model.stochtilfz=tilfz(model.stochfrows,model.stochfzvars);
model.pretilfz=tilfz(model.prefrows,model.prefzvars);

% calculate some big matrices now

if order>=4
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
%     save('tempfile')
%     
%     
%     
%     n_stochfvars=length(stochfvars);
%     clearvars -except n_stochfvars
% 
%     if n_stochfvars>0
%         [ tempM ] = chainsM( n_stochfvars,4 );
%     else
%         tempM=cell(4,1);
%     end
%     load('tempfile')
%     model.stochfvars_chain4c_M2=tempM{2};
%     model.stochfvars_chain4c_M3=tempM{3};
%     model.stochfvars_chain4c_M4=tempM{4};
%     clear tempM
    delete('tempfile.mat')
end

if order>=5
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
%     save('tempfile')
%     
%     
%     
%     n_stochfvars=length(stochfvars);
%     clearvars -except n_stochfvars
% 
%     if n_stochfvars>0
%         [ tempM ] = chainsM( n_stochfvars,5 );
%     else
%         tempM=cell(6,1);
%     end
%     load('tempfile')
%     model.stochfvars_chain5c_M1=tempM{1};
%     model.stochfvars_chain5c_M2=tempM{2};
%     model.stochfvars_chain5c_M3=tempM{3};
%     model.stochfvars_chain5c_M4=tempM{4};
%     model.stochfvars_chain5c_M5=tempM{5};
%     model.stochfvars_chain5c_M6=tempM{6};
% 
%     clear tempM
    delete('tempfile.mat')
end

% totindi=5+8;
% totindi=totindi+(model.pre_n-1)*(1+3+4+6)+1+4+6+9;
% totindi=totindi+3;
% totindi=totindi+4;
% totindi=totindi+2;
% totindi=totindi+2;
% totindi=totindi+1;
% totindi=totindi+9;
% totindi=totindi+8;
% totindi=totindi+8+(model.stoch_n-1)*(1+3+4+6)+1+4+6+9;
% totindi=totindi+2;
% totindi=totindi+9;
% totindi=totindi+7;
% 
% model.totindi=totindi;

rehash;


cd(currentFolder)