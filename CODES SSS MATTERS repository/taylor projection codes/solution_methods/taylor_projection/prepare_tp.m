function model=prepare_tp(f,Phi,yp,y,xp,x,shocks,varargin)
% function model=prepare_tp(f,Phi,yp,y,xp,x,shocks,symparams,order)

% This file prepares codes and data that are used later to compute the
% nonlinear system of Taylor projection and the Jacobian.
% Input args:
%   f = symbolic vector of the model conditions (Ef=0)
%   Phi = symbolic vector of the law of motion of state varialbes.
%   yp,y,xp,x = symbolic vectors of the control and state varaibles for
%   current (y,x) and future (yp,xp) periods.
%   shocks = symbolic shocks.
%   symparams = symbolic vector of the model parameters.
%   order = order of the approximiating Polynomials. Can be 1, 2, 3 or 4.
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
% changed on April 19, 2018

disp('Taylor Projection:')

n_regimes=1;
chi=[];
chip=[];
auxfuns=[];
auxvars=[];
logicalparams=[];
logicalvars=[];

nargs=9;

if nargin==nargs
    symparams=varargin{1};
    order=varargin{2};
elseif nargin==nargs+2
    symparams=varargin{1};
    order=varargin{2};
    auxfuns=varargin{3};
    auxvars=varargin{4};
elseif nargin==nargs+4
    symparams=varargin{1};
    order=varargin{2};
    auxfuns=varargin{3};
    auxvars=varargin{4};
    logicalparams=varargin{5};
    logicalvars=varargin{6};
elseif nargin==nargs+3 
    chip=varargin{1};
    chi=varargin{2};
    symparams=varargin{3};
    order=varargin{4};
    n_regimes=varargin{5};
elseif nargin==nargs+5
    chip=varargin{1};
    chi=varargin{2};
    symparams=varargin{3};
    order=varargin{4};
    n_regimes=varargin{5};
    auxfuns=varargin{6};
    auxvars=varargin{7};
elseif nargin==nargs+7
    chip=varargin{1};
    chi=varargin{2};
    symparams=varargin{3};
    order=varargin{4};
    n_regimes=varargin{5};
    auxfuns=varargin{6};
    auxvars=varargin{7};
    logicalparams=varargin{8};
    logicalvars=varargin{9};
else
    error('wrong number of input arguments')
end

f=f(:);
Phi=Phi(:);
yp=yp(:);
y=y(:);
xp=xp(:);
x=x(:);
shocks=shocks(:);
chip=chip(:);
chi=chi(:);
symparams=symparams(:);
auxfuns=auxfuns(:);
auxvars=auxvars(:);
logicalparams=logicalparams(:);
logicalvars=logicalvars(:);

if length(yp)~=length(y)
    error('current and future control variables must have the same number of elements')
end
if length(xp)~=length(x)
    error('current and future state variables must have the same number of elements')
end
if length(chip)~=length(chi)
    error('current and future Markov-switching parameters must have the same number of elements')
end
if length(auxvars)~=length(auxfuns)
    error('auxiliary functions and variables must have the same number of elements')
end
if length(logicalparams)~=length(logicalvars)
    error('logical parameters and their arguments must have the same number of elements')
end



symparams_ms=chi;
symparams_msp=chip;

fname='ms';
nocheck=0;
uname='';

Phi=sym(Phi);
if isempty(shocks)
    shocks=sym('unused_shock__');
    eta=sym(zeros(length(x),1));
else
    eta=sym(zeros(length(x),length(shocks)));
end

model=prepare_tpMS3(f,Phi,yp,y,xp,x,shocks,symparams_msp,symparams_ms,symparams,eta,...
    n_regimes,order,auxfuns,auxvars,fname,nocheck,uname,logicalparams,logicalvars);
model.pde=0;

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

end