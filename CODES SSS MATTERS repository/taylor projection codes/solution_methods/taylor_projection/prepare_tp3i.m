function model=prepare_tp3i(f,Phi,yp,y,xp,x,shocks,varargin)
% function model=prepare_tp(f,Phi,yp,y,xp,x,symparams,eta,order,varargin)

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
% © Copyright, Oren Levintal, March 19, 2017.

disp('Taylor Projection:')

n_regimes=1;
chi=[];
chip=[];
auxfuns=[];
auxvars=[];
logicalparams=[];
logicalvars=[];

% if isempty(shocks)
%     error('shocks are empty')
% end
nargs=10;

if nargin==nargs
    symparams=varargin{1};
    symparamsi=varargin{2};
    order=varargin{3};
elseif nargin==nargs+2
    symparams=varargin{1};
    symparamsi=varargin{2};
    order=varargin{3};
    auxfuns=varargin{4};
    auxvars=varargin{5};
elseif nargin==nargs+4
    symparams=varargin{1};
    symparamsi=varargin{2};
    order=varargin{3};
    auxfuns=varargin{4};
    auxvars=varargin{5};
    logicalparams=varargin{6};
    logicalvars=varargin{7};
elseif nargin==nargs+3 
    chip=varargin{1};
    chi=varargin{2};
    symparams=varargin{3};
    symparamsi=varargin{4};
    order=varargin{5};
    n_regimes=varargin{6};
elseif nargin==nargs+5
    chip=varargin{1};
    chi=varargin{2};
    symparams=varargin{3};
    symparamsi=varargin{4};
    order=varargin{5};
    n_regimes=varargin{6};
    auxfuns=varargin{7};
    auxvars=varargin{8};
elseif nargin==nargs+7
    chip=varargin{1};
    chi=varargin{2};
    symparams=varargin{3};
    symparamsi=varargin{4};
    order=varargin{5};
    n_regimes=varargin{6};
    auxfuns=varargin{7};
    auxvars=varargin{8};
    logicalparams=varargin{9};
    logicalvars=varargin{10};
else
    error('wrong number of input arguments')
end

n_paramsi=numel(symparamsi);
symparams=[symparams(:);symparamsi(:)];

symparams_ms=chi;
symparams_msp=chip;

fname='ms';
nocheck=1;
uname='';

Phi=sym(Phi);
if isempty(shocks)
    shocks=sym('unused_shock');
    eta=sym(zeros(length(x),1));
else
    eta=sym(zeros(length(x),length(shocks)));
end



model=prepare_tpMS3(f,Phi,yp,y,xp,x,shocks,symparams_msp,symparams_ms,symparams,eta,...
    n_regimes,order,auxfuns,auxvars,fname,nocheck,uname,logicalparams,logicalvars);

fprintf('writing residual function to file...')
home_folder=pwd;
if ~exist([pwd '\files'], 'dir')
   mkdir('files');
end
cd('files')
if isempty(chi)
    make_residual3(f,Phi,yp,y,xp,x,shocks,symparams,order,auxfuns,auxvars,model,logicalparams,logicalvars);
else
    make_residualMS3(f,Phi,yp,y,xp,x,shocks,symparams_msp,symparams_ms,symparams,order,auxfuns,auxvars,model,logicalparams,logicalvars);
end
cd(home_folder);
fprintf('done\n')

model.hybrid=3;

model.n_paramsi=n_paramsi;

model.pde=0;

end