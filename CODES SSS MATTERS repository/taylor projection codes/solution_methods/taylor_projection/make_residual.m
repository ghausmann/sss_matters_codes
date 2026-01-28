function make_residual(f,Phi,yp,y,xp,x,symparams,eta,order,varargin)
% The function creates a Matlab file that computes the residual function. 

if nargin==9
    ufun=[];
    u=[];
    type='normal';
    logicalparams=[];
    logicalvars=[];
elseif nargin==12
    ufun=varargin{1};
    u=varargin{2};
    uu0=varargin{3};
    type='normal';
    logicalparams=[];
    logicalvars=[];
elseif nargin==13
    ufun=varargin{1};
    u=varargin{2};
    uu0=varargin{3};
    type=varargin{4};
    logicalparams=[];
    logicalvars=[];
elseif nargin==14
    ufun=varargin{1};
    u=varargin{2};
    uu0=varargin{3};
    type='normal';
    logicalparams=varargin{4};
    logicalvars=varargin{5};
elseif nargin==15
    ufun=varargin{1};
    u=varargin{2};
    uu0=varargin{3};
    type=varargin{4};
    logicalparams=varargin{5};
    logicalvars=varargin{6};
else
    error('wrong number of input arguments')
end

% search for reserved names

allnames=[yp(:);y(:);xp(:);x(:);symparams(:);u(:)];
reserved_names={'coeffs','x','params','c0','nodes','weights','n_y','n_x1',...
    'n_f','X','n_b','endogenous_vars','y','expected_endogenous_states',...
    'expected_exogenous_states','h_fun','eta_matrix','Resid','i','xp','Xp','yp','f'};

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

% Lenght of vectors
n_f=length(f); n_x=length(x); n_y=length(y); n_x2=size(Phi,1); n_x1=n_x-n_x2;
n_e=size(eta,2); n_u=length(u); 
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
    A{l}='function Resid=residual(coeffs,x,params,c0,nodes,weights)'; l=l+1;
elseif strcmp(type,'MS')
    A=cell(1100+length(symparams)+2*n_x+2*n_y+n_subs*n_u+n_f,1);
    l=1;
    A{l}='function RMS=residual(coeffs,x,params,msparams,transition,c0,nodes,weights)'; l=l+1;

    A{l}='% The function computes the residual function of a model with Markov-switching parameters'; l=l+1;

    A{l}=''; l=l+1;
    A{l}='n_regimes=size(transition,1);'; l=l+1;
    A{l}=['n_f=' num2str(length(f)) ';']; l=l+1;

    A{l}=''; l=l+1;
    A{l}='% compute residual function for all current regimes'; l=l+1;
    A{l}='RMS=zeros(n_f,n_regimes);'; l=l+1;

    A{l}=''; l=l+1;
    A{l}='for i=1:n_regimes'; l=l+1;
    A{l}='   for j=1:n_regimes'; l=l+1;
    A{l}=''; l=l+1;
    A{l}='      % all model parameters given current regime i and future regime j'; l=l+1;
    A{l}='      all_params=[params(:);msparams(:,i);msparams(:,j)];'; l=l+1;

    A{l}=''; l=l+1;
    A{l}='      % model conditions given current regime i and future regime j'; l=l+1;
    A{l}='      R=residual_sub(coeffs(:,i),coeffs(:,j),x,all_params,c0,nodes,weights);'; l=l+1;
    A{l}=''; l=l+1;
    A{l}='      % weight by transition probabilities'; l=l+1;
    A{l}='      RMS(:,i)=RMS(:,i)+R*transition(i,j);'; l=l+1;
    A{l}='   end'; l=l+1;
    A{l}='end'; l=l+1;
    A{l}=''; l=l+1;
    A{l}='end'; l=l+1;
    A{l}=''; l=l+1;
    A{l}=''; l=l+1;
    A{l}='function Resid=residual_sub(coeffs,coeffsp,x,params,c0,nodes,weights)'; l=l+1;
end

A{l}=''; l=l+1;
A{l}='% parameters'; l=l+1;

for i=1:length(symparams)
    A{l}=[char(symparams(i)) '=params(' num2str(i) ');']; l=l+1;
end

A{l}=''; l=l+1;
A{l}='% number of control variables'; l=l+1;
A{l}=['n_y=' num2str(n_y) ';']; l=l+1;

A{l}=''; l=l+1;
A{l}='% number of endogenous state variables'; l=l+1;
A{l}=['n_x1=' num2str(n_x1) ';']; l=l+1;

A{l}=''; l=l+1;
A{l}='% number of all endogenous variables (controls and states)'; l=l+1;
A{l}=['n_f=n_y+n_x1;']; l=l+1;

A{l}=''; l=l+1;
A{l}='% state variables'; l=l+1;

for i=1:n_x
    A{l}=[char(x(i)) '=x(' num2str(i) ');'];l=l+1;
end

A{l}=''; l=l+1;
A{l}='% basis function at the current state x'; l=l+1;

A{l}='X=make_basis(x,c0);'; l=l+1;

A{l}=''; l=l+1;
A{l}='% size of basis function'; l=l+1;

A{l}='n_b=size(X,1);'; l=l+1;



A{l}='coeffs=reshape(coeffs,n_f,n_b);'; l=l+1;


A{l}='endogenous_vars=coeffs*X;'; l=l+1;

A{l}=''; l=l+1;
A{l}='% control variables y'; l=l+1;

A{l}=['y=endogenous_vars(1:n_y);']; l=l+1;

for i=1:n_y
    A{l}=[char(y(i)) '=y(' num2str(i) ');'];l=l+1;
end

A{l}=''; l=l+1;
A{l}='% expected value of future state variables'; l=l+1;

A{l}=['expected_endogenous_states=endogenous_vars(n_y+1:end);']; l=l+1;

A{l}=['expected_exogenous_states=zeros(' num2str(n_x2) ',1);']; l=l+1;

for i=1:n_x2
    A{l}=['expected_exogenous_states(' num2str(i) ')=' char(Phi(i)) ';'];l=l+1;
end

A{l}=['h_fun=[expected_endogenous_states;expected_exogenous_states];']; l=l+1;

A{l}=''; l=l+1;
A{l}='% eta matrix'; l=l+1;

A{l}=['eta_matrix=zeros(' num2str(n_x) ',' num2str(n_e) ');']; l=l+1;

for i=1:n_x
    for j=1:n_e
        if eta(i,j)~=0
            A{l}=['eta_matrix(' num2str(i) ',' num2str(j) ')=' char(eta(i,j)) ';'];l=l+1;
        end
    end
end



A{l}=''; l=l+1;
A{l}='% compute the residual function'; l=l+1;

A{l}='Resid=zeros(n_f,1);'; l=l+1;

A{l}=''; l=l+1;
A{l}='% use quadrature nodes and weights to approximate expectations'; l=l+1;

A{l}='for i=1:length(weights)'; l=l+1;

A{l}=''; l=l+1;
A{l}='   % future state variables for node i'; l=l+1;

A{l}='   xp=h_fun+eta_matrix*nodes(:,i);'; l=l+1;

for i=1:n_x
    A{l}=['   ' char(xp(i)) '=xp(' num2str(i) ');'];l=l+1;
end

A{l}=''; l=l+1;
A{l}='   % basis function at the future state xp'; l=l+1;

A{l}='   Xp=make_basis(xp,c0);'; l=l+1;

if strcmp(type,'normal')
    A{l}='   yp=coeffs(1:n_y,:)*Xp;'; l=l+1;
elseif strcmp(type,'MS')
    A{l}='   coeffsp=reshape(coeffsp,n_f,n_b);'; l=l+1;
    A{l}='   yp=coeffsp(1:n_y,:)*Xp;'; l=l+1;
end

for i=1:n_y
    A{l}=['   ' char(yp(i)) '=yp(' num2str(i) ');'];l=l+1;
end

if ~isempty(u)
    A{l}=''; l=l+1;
    A{l}='   % auxiliary variables'; l=l+1;

    % v=[yp(:);y(:);xp(:);x(:)];


    % subsu=ufun;
    % newsubsu=subs(subsu,u,ufun);
    % 
    % uind=1:n_u;
    % while ~isequal(subsu,newsubsu)
    %     for i=uind(logical(newsubsu==subsu))
    %         A{l}=['   ' char(u(i)) '=' char(ufun(i)) ';'];l=l+1;
    %     end
    %     uind=uind(logical(newsubsu~=subsu));
    %     subsu=newsubsu(logical(newsubsu~=subsu));
    %     newsubsu=subs(subsu,u,ufun);
    % end
    % for i=uind
    %     A{l}=['   ' char(u(i)) '=' char(ufun(i)) ';'];l=l+1;
    % end

    % version 2

    % if nargin~=8
%         uu0=double(logical(jacobian(ufun(:),u)~=0));
    % end

    
    uu=uu0;

    n_u=length(u);
    uind=1:n_u;
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
end

%%%%%%%%%%

if ~isempty(logicalparams)
    A{l}=''; l=l+1;
    for i=1:length(logicalparams)
        A{l}=['   ' char(logicalparams(i)) '=double(logical(' char(logicalvars(i)) '>=0));'];l=l+1;
    end
end

A{l}=''; l=l+1;
A{l}=['   % compute function f'];l=l+1;

A{l}=['   f=zeros(n_f,1);'];l=l+1;
for i=1:n_f
    A{l}=['   f(' num2str(i) ')=' char(f(i)) ';'];l=l+1;
end

A{l}=''; l=l+1;
A{l}=['   % multiply by weight i and add to Resid'];l=l+1;

A{l}=['   Resid=Resid+f*weights(i);'];l=l+1;

A{l}='end'; l=l+1;
A{l}=''; l=l+1;
A{l}='end'; l=l+1;


% Write cell A into txt

fid = fopen(['residual.m'], 'w');
for i = 1:l-1
    fprintf(fid,'%s\n', A{i});
end

fclose(fid);

rehash;

