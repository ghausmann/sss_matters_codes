function make_residual3(f,Phi,yp,y,xp,x,shocks,symparams,order,varargin)
% The function creates a Matlab file that computes the residual function. 

nargs=0;

if nargin==nargs+9
    ufun=[];
    u=[];
    type='normal';
    logicalparams=[];
    logicalvars=[];
elseif nargin==nargs+12
    ufun=varargin{1};
    u=varargin{2};
    model=varargin{3};
    type='normal';
    logicalparams=[];
    logicalvars=[];
elseif nargin==nargs+13
    ufun=varargin{1};
    u=varargin{2};
    model=varargin{3};
    type=varargin{4};
    logicalparams=[];
    logicalvars=[];
elseif nargin==nargs+14
    ufun=varargin{1};
    u=varargin{2};
    model=varargin{3};
    type='normal';
    logicalparams=varargin{4};
    logicalvars=varargin{5};
elseif nargin==nargs+15
    ufun=varargin{1};
    u=varargin{2};
    model=varargin{3};
    type=varargin{4};
    logicalparams=varargin{5};
    logicalvars=varargin{6};
else
    error('wrong number of input arguments')
end

uu0=model.piu0;
tilPhiu=model.tilPhiu*model.uu; % this is a matrix for all auxiliary variables that affect Phi (through all substitutions)



% search for reserved names

allnames=[yp(:);y(:);xp(:);x(:);symparams(:);u(:);shocks(:)];
reserved_names={'coeffs','coeffsp','state','future_state','control','future_control','regime','future_regime','params','msparams','regimes','n_regimes',...
    'shock','center','nodes','n_nodes','weights','Basis','n_nodes','endogenous_vars','transition',...
    'R_fun','g_fun','Phi_fun','aux_fun','i','j','f_fun',...
    'R_funMS','g_funMS','Phi_funMS','aux_funMS'};

reserved='';


for i=1:length(reserved_names)
    if nnz(logical(allnames==sym(reserved_names{i})))~=0
        reserved=[reserved, ' \n ', reserved_names{i}];
    end
end

if ~isempty(reserved)
    error_message=['The following names are reserved and should be changed:', reserved];
    error('ErrorTests:convertTest',error_message)
end

% Length of vectors
n_f=length(f); n_x=length(x); n_y=length(y); n_x2=size(Phi,1); n_x1=n_x-n_x2;
n_e=length(shocks); n_u=length(u); 
% n_subs=find_n_new(ufun,u);

n_subs=1;
if ~isempty(u)
    uu=sparse(uu0);
    while nnz(uu)>0
        n_subs=n_subs+1;
        uu=uu*uu0;
    end
end
order=order(1);
n_b=nchoosek(n_x+order,order);

% make a basis function

A=cell(1000+n_b+n_x,1);
l=1;

A{l}='function X=make_basis(x,c0)'; l=l+1;

A{l}='% The function computes the basis function at x'; l=l+1;

A{l}=''; l=l+1;
A{l}='% compute x-c0'; l=l+1;

x_c=sym(zeros(n_x,1));
for i=1:n_x
    A{l}=['x' num2str(i) '=x(' num2str(i) ')-c0(' num2str(i) ');']; l=l+1;
    x_c(i)=sym(['x' num2str(i)]);
end

A{l}=''; l=l+1;
A{l}='% create the basis function X'; l=l+1;


kronx=x_c;
symX=sym(1);
if order>=1
    symX=[sym(1);kronx];
end
for i=2:order
    kronx=kron(x_c,kronx);
    [~,W]=create_UW(n_x,i);
    symX=[symX;W*kronx];
end

A{l}=['X=zeros(' num2str(length(symX)) ',1);']; l=l+1;

for i=1:length(symX)
    A{l}=['X(' num2str(i) ')=' char(symX(i)) ';'];
    l=l+1;
end

% Write cell A into txt

fid = fopen(['make_basis.m'], 'w');
for i = 1:l-1
    fprintf(fid,'%s\n', A{i});
end

fclose(fid);

%%%%%%%%%%%%%%%%%%%%

if strcmp(type,'normal')
    A=cell(1000+length(symparams)+2*n_x+2*n_y+n_subs*n_u+n_f,1);
    l=1;
    A{l}='function [R_fun,g_fun,Phi_fun,aux_fun]=residual(coeffs,state,params,center,nodes,weights)'; l=l+1;
    A{l}='% The function computes the residual function for a model with iid shocks'; l=l+1;
elseif strcmp(type,'MS')
    A=cell(1100+length(symparams)+2*n_x+2*n_y+n_subs*n_u+n_f,1);
    l=1;
    A{l}='function [R_funMS,g_funMS,Phi_funMS,aux_funMS]=residual(coeffs,state,params,msparams,transition,center,nodes,weights,varargin)'; l=l+1;

    A{l}='% The function computes the residual function for a model with iid shocks and Markov shocks'; l=l+1;

    A{l}=''; l=l+1;
    
    A{l}='if nargin==8'; l=l+1;
    A{l}=['   regimes=1:' num2str(model.n_regimes) ';']; l=l+1;
    A{l}='else'; l=l+1;
    A{l}=['   regimes=varargin{1};']; l=l+1;
    A{l}='end'; l=l+1;
    A{l}='n_regimes=length(regimes);'; l=l+1;    
    A{l}=''; l=l+1;    
%     A{l}=['n_f=' num2str(length(f)) ';']; l=l+1;
%     A{l}=['n_y=' num2str(length(y)) ';']; l=l+1;
% %     A{l}=['n_x1=' num2str(n_x1) ';']; l=l+1;
%     A{l}=['n_x=' num2str(n_x2) ';']; l=l+1;
    A{l}=['n_nodes=length(weights);']; l=l+1;
%     A{l}='n_regimes=size(transition,1);'; l=l+1;

    A{l}=''; l=l+1;
    A{l}=['coeffs=reshape(coeffs,' num2str(model.n_y*model.n_b) ',' num2str(model.n_regimes) ');']; l=l+1;
    A{l}=['msparams=reshape(msparams,' num2str(model.n_msparams) ',' num2str(model.n_regimes) ');']; l=l+1;
    
    A{l}=''; l=l+1;
    A{l}='% preallocate variables'; l=l+1;
    A{l}=['R_funMS=zeros(' num2str(model.n_f) ',n_regimes);']; l=l+1;
    A{l}=['g_funMS=zeros(' num2str(model.n_y) ',n_regimes);']; l=l+1;
%     A{l}='h_funMS=zeros(n_x1,n_regimes);'; l=l+1;
    A{l}=['Phi_funMS=zeros(' num2str(model.n_x) ',n_nodes,n_regimes,' num2str(model.n_regimes) ');']; l=l+1;
    A{l}=['aux_funMS=zeros(' num2str(n_u) ',n_nodes,n_regimes,' num2str(model.n_regimes) ');']; l=l+1;

    A{l}=''; l=l+1;
    A{l}='loc=0;'; l=l+1;
    A{l}=['for i=regimes']; l=l+1;
    A{l}=['   loc=loc+1;']; l=l+1;
    A{l}=['   for j=1:' num2str(model.n_regimes)]; l=l+1;
    A{l}=''; l=l+1;
    A{l}='      % all model parameters given current regime i and future regime j'; l=l+1;
    A{l}='      all_params=[params(:);msparams(:,i);msparams(:,j)];'; l=l+1;

    A{l}=''; l=l+1;
    A{l}='      % model conditions given current regime i and future regime j'; l=l+1;
    A{l}='      [R_fun,g_fun,Phi_fun,aux_fun]=residual_sub(coeffs,state,i,j,all_params,center,nodes,weights);'; l=l+1;
    A{l}=''; l=l+1;
    A{l}='      % weight by transition probabilities'; l=l+1;
    A{l}='      R_funMS(:,loc)=R_funMS(:,loc)+R_fun*transition(i,j);'; l=l+1;
    A{l}='      Phi_funMS(:,:,loc,j)=Phi_fun;'; l=l+1;
    A{l}='      aux_funMS(:,:,loc,j)=aux_fun;'; l=l+1;
    A{l}='   end'; l=l+1;
    A{l}='   g_funMS(:,loc)=g_fun;'; l=l+1;
%     A{l}='   h_funMS(:,i)=reshape(h_fun,n_x1,1);'; l=l+1;
    
    A{l}='end'; l=l+1;
    A{l}=''; l=l+1;
    A{l}='end'; l=l+1;
    A{l}=''; l=l+1;
    A{l}=''; l=l+1;
    A{l}='function [R_fun,g_fun,Phi_fun,aux_fun]=residual_sub(coeffs,state,regime,future_regime,params,center,nodes,weights)'; l=l+1;
end

l1=l;
A{l}=''; l=l+1;
A{l}='% parameters'; l=l+1;

for i=1:length(symparams)
    A{l}=[char(symparams(i)) '=params(' num2str(i) ');']; l=l+1;
end
l2=l;

% A{l}=''; l=l+1;
% A{l}='% number of control variables'; l=l+1;
% A{l}=['n_y=' num2str(n_y) ';']; l=l+1;

% A{l}=''; l=l+1;
% A{l}='% number of predetermined state variables'; l=l+1;
% A{l}=['n_x1=' num2str(n_x1) ';']; l=l+1;
% 
% A{l}=''; l=l+1;
% A{l}='% number of endogenous variables'; l=l+1;
% A{l}=['n_endog=n_y+n_x1;']; l=l+1;

A{l}=''; l=l+1;
l3=l;

A{l}='% state variables'; l=l+1;

for i=1:n_x
    A{l}=[char(x(i)) '=state(' num2str(i) ');'];l=l+1;
end
l4=l;

% A{l}='% basis function at the current state'; l=l+1;
% 
% A{l}='Basis=make_basis(state,c0);'; l=l+1;
% 
% A{l}=''; l=l+1;
% A{l}='% size of basis function'; l=l+1;
% 
% A{l}='n_b=size(Basis,1);'; l=l+1;


% A{l}=''; l=l+1;

% A{l}=['coeffs=reshape(coeffs,' num2str(model.n_y) ',' num2str(model.n_b) ');']; l=l+1;
% if strcmp(type,'MS')
%     A{l}='coeffsp=reshape(coeffsp,n_y,n_b);'; l=l+1;
% end

% A{l}='endogenous_vars=coeffs*X;'; l=l+1;

A{l}=''; l=l+1;
A{l}='% the function g(x)'; l=l+1;
if strcmp(type,'normal')
    A{l}=['g_fun=evalg(state,coeffs,center);']; l=l+1;
else
    A{l}=['g_fun=evalg(state,regime,coeffs,center);']; l=l+1;
end

l5=l;
A{l}=''; l=l+1;
A{l}='% control variables'; l=l+1;


for i=1:n_y
    A{l}=[char(y(i)) '=g_fun(' num2str(i) ');'];l=l+1;
end
l6=l;
%%%%%



%%%%%

A{l}=''; l=l+1;
A{l}=['if isempty(weights) || isempty(nodes)']; l=l+1;
A{l}=['   weights=1;']; l=l+1;
A{l}=['   nodes=0;']; l=l+1;
A{l}=['end']; l=l+1;

A{l}=''; l=l+1;
A{l}=['n_nodes=length(weights);']; l=l+1;
% A{l}=['n_f=' num2str(length(f)) ';']; l=l+1;



% A{l}=['h_fun=[expected_endogenous_states;expected_exogenous_states];']; l=l+1;
% A{l}=['h_fun=future_endogenous_states;']; l=l+1;

A{l}=''; l=l+1;
A{l}='% preallocate variables'; l=l+1;
A{l}=['f_fun=zeros(' num2str(model.n_f) ',1);'];l=l+1;
A{l}=['R_fun=zeros(' num2str(model.n_f) ',1);']; l=l+1;
A{l}=['Phi_fun=zeros(' num2str(n_x2) ',n_nodes);']; l=l+1;
A{l}=['aux_fun=zeros(' num2str(n_u) ',n_nodes);']; l=l+1;


A{l}=''; l=l+1;
A{l}='% approximate expectations by quadrature nodes and weights'; l=l+1;

A{l}='for i=1:n_nodes'; l=l+1;

A{l}=''; l=l+1;
A{l}='   % realization of shocks at node i'; l=l+1;

for j=1:length(shocks)
    A{l}=['   ' char(shocks(j)) '=nodes(' num2str(j) ',i);'];l=l+1;
end

if isempty(u)
    A{l}=''; l=l+1;
    for i=1:n_x2
        A{l}=['   Phi_fun(' num2str(i) ',i)=' char(Phi(i)) ';'];l=l+1;
    end
else
    l7=l;
    A{l}=''; l=l+1;
    
    
    n_u=length(u);
    uind=1:n_u;
    
    uu=uu0;
    newuind=uind;
    Phiuind=find(logical(sum(tilPhiu,1)~=0));
    while nnz(uu(Phiuind,:))~=0
        
        for i=Phiuind(logical(sum(uu(Phiuind,:),2)==0))
            if l==l7+1
                A{l}='   % auxiliary variables'; l=l+1;
            end
            A{l}=['   ' char(u(i)) '=' char(ufun(i)) ';'];l=l+1;
            newuind(i)=0;
        end
        Phiuind=Phiuind(logical(sum(uu(Phiuind,:),2)~=0));
        uu=uu*uu0;
    end
    for i=Phiuind
        A{l}=['   ' char(u(i)) '=' char(ufun(i)) ';'];l=l+1;
        newuind(i)=0;
    end
    l8=l;
    A{l}=''; l=l+1;
    A{l}='   % function Phi'; l=l+1;

    for i=1:n_x2
        A{l}=['   Phi_fun(' num2str(i) ',i)=' char(Phi(i)) ';'];l=l+1;
    end
    
%     A{l}=''; l=l+1;
%     A{l}='   % store auxiliary variables'; l=l+1;
%     for i=1:length(u)
%         A{l}=['   auxvars(' num2str(i) ',i)=' char(u(i)) ';'];l=l+1;
%     end
end

A{l}=''; l=l+1;
A{l}='   % future state variables'; l=l+1;

A{l}='   future_state=Phi_fun(:,i);'; l=l+1;

% A{l}=''; l=l+1;
for i=1:n_x
    A{l}=['   ' char(xp(i)) '=future_state(' num2str(i) ');'];l=l+1;
end

% A{l}=''; l=l+1;
% A{l}='   % basis function at the future state xp'; l=l+1;
% 
% A{l}='   Xp=make_basis(xp,c0);'; l=l+1;

A{l}=''; l=l+1;
A{l}='   % control variables at the future state'; l=l+1;

if strcmp(type,'normal')
    A{l}='   future_control=evalg(future_state,coeffs,center);'; l=l+1;
elseif strcmp(type,'MS')
    A{l}='   future_control=evalg(future_state,future_regime,coeffs,center);'; l=l+1;
end

for i=1:n_y
    A{l}=['   ' char(yp(i)) '=future_control(' num2str(i) ');'];l=l+1;
end

if ~isempty(u) && nnz(newuind)>0
    A{l}=''; l=l+1;
    A{l}='   % auxiliary variables'; l=l+1;

    
    uu=uu0;

    uind=nonzeros(newuind);
    uind=uind(:)';
    uu=uu(uind,:);
    while nnz(uu)~=0
        for i=uind(logical(sum(uu,2)==0))
            A{l}=['   ' char(u(i)) '=' char(ufun(i)) ';'];l=l+1;
        end
        uind=uind(logical(sum(uu,2)~=0));
        uu=uu(logical(sum(uu,2)~=0),:);
        uu=uu*uu0;
    end
    for i=uind
        A{l}=['   ' char(u(i)) '=' char(ufun(i)) ';'];l=l+1;
    end
    
    A{l}=''; l=l+1;
    A{l}='   % store auxiliary variables'; l=l+1;
    for i=1:n_u
        A{l}=['   aux_fun(' num2str(i) ',i)=' char(u(i)) ';'];l=l+1;
    end
    
end

%%%%%%%%%%

if ~isempty(logicalparams)
    A{l}=''; l=l+1;
    A{l}=['   % logical parameters'];l=l+1;
    for i=1:length(logicalparams)
        A{l}=['   ' char(logicalparams(i)) '=double(logical(' char(logicalvars(i)) '>0));'];l=l+1;
    end
end

A{l}=''; l=l+1;
A{l}=['   % compute function f'];l=l+1;

for i=1:n_f
    A{l}=['   f_fun(' num2str(i) ')=' char(f(i)) ';'];l=l+1;
end

A{l}=''; l=l+1;
A{l}=['   % multiply by weight i and add to R_fun'];l=l+1;

A{l}=['   R_fun=R_fun+f_fun*weights(i);'];l=l+1;

A{l}='end'; l=l+1;
A{l}=''; l=l+1;
A{l}='end'; l=l+1;


% Write cell A into txt

fid = fopen(['residual.m'], 'w');
for i = 1:l-1
    fprintf(fid,'%s\n', A{i});
end

fclose(fid);


% Write evalg

fid = fopen(['evalg.m'], 'w');

if strcmp(type,'normal')
    fprintf(fid,'%s\n', 'function g_fun=evalg(x,coeffs,c0)');
    fprintf(fid,'%s\n', '% The function computes the function g(x), where x is the state vector');
else
    fprintf(fid,'%s\n', 'function g_fun=evalg(x,s,coeffs,c0)');
    fprintf(fid,'%s\n', '% The function computes the function g(x,s), where x is the state vector and s is the regime');
end

fprintf(fid,'%s\n', '');

fprintf(fid,'%s\n', ['Basis=make_basis(x,c0);']);
fprintf(fid,'%s\n', '');
if strcmp(type,'MS')
    fprintf(fid,'%s\n', ['coeffs=reshape(coeffs,' num2str(model.n_y*model.n_b) ',' num2str(model.n_regimes) ');']);
    fprintf(fid,'%s\n', ['coeffs=coeffs(:,s);']);
end
fprintf(fid,'%s\n', ['coeffs=reshape(coeffs,' num2str(model.n_y) ',' num2str(model.n_b) ');']);
fprintf(fid,'%s\n', '');
fprintf(fid,'%s\n', ['g_fun=coeffs*Basis;']);
fprintf(fid,'%s\n', '');

fprintf(fid,'%s\n', ['end']);


fclose(fid);


% Write evalPhi

fid = fopen(['evalPhi.m'], 'w');

if strcmp(type,'normal')
    fprintf(fid,'%s\n', 'function Phi_fun=evalPhi(state,control,shock,params)');
    fprintf(fid,'%s\n', '% The function computes the function Phi(x,y,epsp), where x is the state vector,');
    fprintf(fid,'%s\n', '% y is the control vector given by y=g(x), and epsp is the realization of the shocks next period.');
else
    fprintf(fid,'%s\n', 'function Phi_fun=evalPhi(state,control,shock,future_regime,regime,params,msparams)');
    fprintf(fid,'%s\n', '% The function computes the function Phi(x,y,epsp,chi(sp),chi(s)), where x is the state vector,');
    fprintf(fid,'%s\n', '% y is the control vector given by y=g(x,s), chi(s) is the value of the Markov shocks at the current regime s,');
    fprintf(fid,'%s\n', '% chi(sp) is the value of the Markov shocks at the future regime sp, and epsp is the realization of the shocks next period.');
end

% parameters
if strcmp(type,'MS')
    fprintf(fid,'%s\n', '');
    temp='params=[params(:);msparams(:,regime);msparams(:,future_regime)];'; 
    fprintf(fid,'%s\n', temp);
end

for i = l1:l2
    fprintf(fid,'%s\n', A{i});
end


% state variables
for i = l3:l4
    fprintf(fid,'%s\n', A{i});
end

% control variables
temp=['% control variables'];
fprintf(fid,'%s\n', temp);

for j=1:model.n_y
    temp=[char(y(j)) '=control(' num2str(j) ');'];
    fprintf(fid,'%s\n', temp);
end
fprintf(fid,'%s\n', '');

temp=['% shocks'];
fprintf(fid,'%s\n', temp);

% shocks
for j=1:length(shocks)
    temp=[char(shocks(j)) '=shock(' num2str(j) ');'];
    fprintf(fid,'%s\n', temp);
end
fprintf(fid,'%s\n', '');

temp=['% function Phi'];
fprintf(fid,'%s\n', temp);

if ~isempty(u)
    for i = l7:l8
        fprintf(fid,'%s\n', A{i}(4:end));
    end
end

temp=['Phi_fun=zeros(' num2str(model.n_x2) ',1);'];
fprintf(fid,'%s\n', temp);

for i=1:n_x2
    temp=['Phi_fun(' num2str(i) ')=' char(Phi(i)) ';'];
    fprintf(fid,'%s\n', temp);
end

fprintf(fid,'%s\n', '');

fprintf(fid,'%s\n', ['end']);


fclose(fid);

rehash;

