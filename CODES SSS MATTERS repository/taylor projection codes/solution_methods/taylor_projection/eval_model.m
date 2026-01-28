function [R0,g,h1,nPhi,nu_vec,model]=eval_model(coeffs,x0,model,params,varargin)
% function [R0,g,h1,nPhi,nu_vec,model]=eval_model(coeffs,x0,model,params,c0,nep,P)

if nargin==7
    type='normal';
    c0=varargin{1};
    nep=varargin{2};
    P=varargin{3};
elseif nargin==9
    type='MS';
    msparams=varargin{1};
    transition=varargin{2};
    c0=varargin{3};
    nep=varargin{4};
    P=varargin{5};
    specific_regime=[];
elseif nargin==10
    type='MS';
    msparams=varargin{1};
    transition=varargin{2};
    c0=varargin{3};
    nep=varargin{4};
    P=varargin{5};
    specific_regime=varargin{6};
else
    error('wrong number of input arguments')
end

if isempty(nep) || isempty(P)
    P=1;
    nep=0;
end

if strcmp(type,'normal')

    coeffs=coeffs(:);
    if length(coeffs)~=model.n_theta
        error('wrong number of polynomial coefficients')
    end
    x0=x0(:);
    if length(x0)~=model.n_x
        error('wrong size of approximation point')
    end

%     params=[params(:);zeros(model.n_logicalparams,1)];

    if length(params)~=model.n_params-model.n_logicalparams
        error('wrong number of model parameters')
    end
    c0=c0(:);
    if length(c0)~=model.n_x
        error('wrong size of polynomial center')
    end
    if size(nep,1)~=model.n_e
        error('wrong number of shocks')
    end
    P=P(:);
    if size(nep,2)~=length(P)
        error('incorrect distribution of shocks')
    end

    msparams=zeros(0,1);
    transition=1;
    [R0,g,h1,nPhi,nu_vec,model]=eval_modelMS(coeffs,x0,model,params,msparams,transition,c0,nep,P);

    g=reshape(g,model.n_y,1);
    h1=reshape(h1,model.n_x1,1);
    nPhi=reshape(nPhi,model.n_x2,[]);
    n_nodes=length(P(:));
    nu_vec=reshape(nu_vec,model.n_u,n_nodes); % n_u,n_nodes

else

    if size(coeffs,1)~=model.n_theta || size(coeffs,2)~=model.n_regimes
        error('wrong number of polynomial coefficients')
    end
    x0=x0(:);
    if length(x0)~=model.n_x
        error('wrong size of approximation point')
    end
    
%     params=[params(:);zeros(model.n_logicalparams,1)];

    if length(params)~=model.n_params_no_ms-model.n_logicalparams
        error('wrong number of model parameters')
    end
    if isempty(msparams)
        msparams=zeros(0,model.n_regimes);
    end
    if size(msparams,1)~=model.n_msparams || size(msparams,2)~=model.n_regimes
        error('wrong size of Markov switching parameters')
    end
    if size(transition,1)~=model.n_regimes || size(transition,2)~=model.n_regimes
        error('wrong size of transition matrix')
    end
    c0=c0(:);
    if length(c0)~=model.n_x
        error('wrong size of polynomial center')
    end
    if size(nep,1)~=model.n_e
        error('wrong number of shocks')
    end
    P=P(:);
    if size(nep,2)~=length(P)
        error('incorrect distribution of shocks')
    end    
    
    if isempty(specific_regime)
        [R0,g,h1,nPhi,nu_vec,model]=eval_modelMS(coeffs,x0,model,params,msparams,transition,c0,nep,P);
    else
        [R0,g,h1,nPhi,nu_vec,model]=eval_modelMS(coeffs,x0,model,params,msparams,transition,c0,nep,P,specific_regime);
    end
    
    
end