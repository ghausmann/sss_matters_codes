function [ model ] = prepare_tpMS2( f,Phi,yp,y,xp,x,shocks,symparams_msp,symparams_ms,symparams,...
    eta,n_regimes,order,auxfuns,auxvars,varargin )


n_f=length(f); n_x=length(x); n_y=length(y); n_x2=length(Phi); n_x1=n_x-n_x2;

fname='ms';
nocheck=0;
uname='';
logicalparams=[];
logicalvars=[];

if ~isempty(varargin)
    fname=varargin{1};
    if ~ischar(fname)
        error('function name must be a string')
    end    
    if length(varargin)>1
        nocheck=varargin{2};
    end
    if length(varargin)>2
        uname=varargin{3};
    end
    if length(varargin)>3
        logicalparams=varargin{4};
        logicalvars=varargin{5};
    end
end

if nocheck==0
    fprintf('preprocessing (this may take a while)...')
end

if n_f~=n_y+n_x1 && nocheck==0
    error('number of model conditions should be equal to the number of endogenous variables')
end

if numel(symparams_msp)~=numel(symparams_ms)
    error('wrong size of future and current Markov switching parameters')
end

new_symparams=[symparams(:);symparams_ms(:);symparams_msp(:)];

model=prepare_tpvec2(f,Phi,yp,y,xp,x,shocks,new_symparams,eta,order,auxfuns,auxvars,fname,uname,logicalparams,logicalvars);
model.n_regimes=n_regimes;

model.n_msparams=numel(symparams_ms);
model.n_params_no_ms=numel(symparams)+length(logicalparams);

model.steady_state=0;

if nocheck==0
    fprintf('done\n')
end

end

