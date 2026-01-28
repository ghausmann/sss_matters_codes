function [R0,g,h1,nPhi,nu_vec,model]=eval_modelMS(coeffs,x0,model,params,msparams,transition,c0,nep,P,varargin)

x0=x0(:);

params=params(:);

c0=c0(:);

P=P(:);

n_nodes=numel(P);

if isempty(varargin)

    msparams=msparams';

    paramsvec=repmat(params(:)',model.n_regimes^2,1); % n_regimes*n_regimes,[]

    [s,sp]=ndgrid([1:model.n_regimes],[1:model.n_regimes]);
    s=s(:);
    sp=sp(:);
    paramsvec=[paramsvec,msparams(s,:),msparams(sp,:)]; %n_s,n_sp,n_params

    paramsvec=[paramsvec,zeros(model.n_regimes^2,model.n_logicalparams)];
    
    coeffs=reshape(full(coeffs),model.n_theta,model.n_regimes);
    coeffsvec=coeffs(:,s)';
    coeffspvec=coeffs(:,sp)';



    x0vec=repmat(x0(:)',model.n_regimes^2,1);
    c0vec=repmat(c0(:)',model.n_regimes^2,1);

    nep=full(nep);
    P=full(P);

    [R0,g,h1,nPhi,nu_vec,model]=eval_model_vec(coeffsvec,x0vec,model,paramsvec,c0vec,nep,P,coeffspvec);

    ny=reshape(g.vals,model.n_regimes,model.n_regimes,[]);
    ny=ny(:,1,:);
    g=reshape(ny,model.n_regimes,model.n_y)';

    nx1=reshape(h1.vals,model.n_regimes,model.n_regimes,[]);
    nx1=nx1(:,1,:);
    h1=reshape(nx1,model.n_regimes,model.n_x1)';

    nPhi=reshape(nPhi,model.n_regimes,model.n_regimes,[]);
    nPhi=nPhi(:,1,:);
    if model.hybrid==0
        nPhi=reshape(nPhi,model.n_regimes,model.n_x2)';
    else
        nPhi=permute(reshape(nPhi,model.n_regimes,model.n_x2,n_nodes),[2,1,3]);
    end
    
    n_nodes=length(P);
    nu_vec=reshape(nu_vec,n_nodes,model.n_regimes,model.n_regimes,model.n_u); % n_nodes,s,sp,n_u
    nu_vec=permute(nu_vec,[4,2,3,1]);
    
    % sum across sp
    R0=reshape(R0,model.n_regimes,model.n_regimes,[]);
    R0=R0.*reshape(repmat(transition,1,size(R0,3)),model.n_regimes,model.n_regimes,size(R0,3));
    R0=reshape(sum(R0,2),model.n_regimes,[])';

    R0=R0(:);

else
    
    regime=varargin{1};
    if numel(regime)>1
        error('only 1 regime is allowed')
    end

    msparams=msparams';

    paramsvec=repmat(params(:)',model.n_regimes,1); % n_regimes*n_regimes,[]

    s=repmat(regime,model.n_regimes,1);
    sp=(1:model.n_regimes)';
    paramsvec=[paramsvec,msparams(s,:),msparams(sp,:)]; %n_s,n_sp,n_params

    paramsvec=[paramsvec,zeros(model.n_regimes,model.n_logicalparams)];
    
    coeffs=reshape(full(coeffs),model.n_theta,model.n_regimes);
    coeffsvec=coeffs(:,s)';
    coeffspvec=coeffs(:,sp)';

    x0vec=repmat(x0(:)',model.n_regimes,1);
    c0vec=repmat(c0(:)',model.n_regimes,1);

    nep=full(nep);
    P=full(P);

    [R0,g,h1,nPhi,nu_vec,model]=eval_model_vec(coeffsvec,x0vec,model,paramsvec,c0vec,nep,P,coeffspvec);

    ny=reshape(g.vals,1,model.n_regimes,[]);
    ny=ny(:,1,:);
    g=squeeze(ny);

    nx1=reshape(h1.vals,1,model.n_regimes,[]);
    nx1=nx1(:,1,:);
    h1=squeeze(nx1);

    if model.hybrid==0
        nPhi=reshape(nPhi,1,model.n_regimes,[]);
        nPhi=nPhi(:,1,:);
        nPhi=squeeze(nPhi);
    else
        nPhi=permute(reshape(nPhi,model.n_regimes,model.n_x2,n_nodes),[2,1,3]);
%         nPhi=nPhi(:,1,:);
        nPhi=reshape(nPhi,model.n_x2,model.n_regimes,n_nodes);
    end
    
    n_nodes=length(P);
    nu_vec=reshape(nu_vec,n_nodes,model.n_regimes,model.n_u); % n_nodes,sp,n_u
    nu_vec=permute(nu_vec,[3,2,1]);
    % sum across sp
    R0=reshape(R0,1,model.n_regimes,[]);
    R0=R0.*reshape(repmat(transition(regime,:),1,size(R0,3)),1,model.n_regimes,size(R0,3));

    R0=reshape(sum(R0,2),1,[]);

    R0=R0(:);
    
end

end