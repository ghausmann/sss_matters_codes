function [ scaledT,scaledJ,Wcols,Wrows,model ] = tpscaleMSi( weighted_coeffs,x0,model,params,paramsi,msparams,transition,c0,nep,P,Wcols,Wrows)
%This function computes the scaled system and its Jacobian. If weights are
%empty, the function computes and returns the weights. Otherwise, the
%supplied weights are used.
%
% © Copyright, Oren Levintal, June 13, 2016.

order=model.order(1);

if isempty(Wcols) && isempty(Wrows) 
    % if weights are not supplied assume no weights.
    coeffs=weighted_coeffs;
elseif ~isempty(Wcols) && ~isempty(Wrows)
    % if all weights are supplied, unscale weighted_coeffs
    coeffs=weighted_coeffs./Wcols;
else
    error('missing weights');
end

if ~isempty(paramsi)
    n_s=size(paramsi,2);
else
    n_s=1;
end

if ~isfield(model,'active_vars')
    coeffs=repmat(coeffs(:)',n_s,1); 
    coeffs=reshape(coeffs,n_s,model.n_theta,model.n_regimes);
else
    coeffs=reshape(coeffs,[],model.n_b);
    coeffs=coeffs(model.active_vars(:),:);
    coeffs=reshape(coeffs,n_s,[]);
end
x0=repmat(x0(:)',n_s,1);
c0=repmat(c0(:)',n_s,1);

params=[repmat(params(:)',n_s,1),paramsi'];

[Ti,Ji,model]=tpMSvec(coeffs,x0,model,params,msparams,transition,c0,nep,P);

if ~isfield(model,'active_vars')
    T=Ti;
    T.vals=Ti.vals(1,:);
    J=Ji;
    J.vals=Ji.vals(1,:);

    for j=2:n_s
        Tj=Ti;
        Tj.vals=Ti.vals(j,:);
        T=vconcat(T,Tj);
        Jj=Ji;
        Jj.vals=Ji.vals(j,:);
        J=vconcat(J,Jj);
    end
else
    %%%%%%%%%%%%%%%
    %%% The coefficients must be of dimensions model.n_vars,n_b, where model.n_vars is
    %%% the number of variables, and n_b size
    %%% of basis function
    %%%%%%%%%%%%%%%
    T=Ti;
    T.vals=Ti.vals(1,:);
    J=Ji;
    J.vals=Ji.vals(1,:);
    countcoeffs=reshape(1:model.n_vars*model.n_b,[],model.n_b);
    J=changecols(J,reshape(countcoeffs(model.active_vars(1,:),:),1,[]),model.n_vars*model.n_b,1);
    for j=2:n_s
        Tj=Ti;
        Tj.vals=Ti.vals(j,:);
        T=vconcat(T,Tj);
        Jj=Ji;
        Jj.vals=Ji.vals(j,:);
        Jj=changecols(Jj,reshape(countcoeffs(model.active_vars(j,:),:),1,[]),model.n_vars*model.n_b,1);
        J=vconcat(J,Jj);
    end
end
T=sptensor2spmat(T);
J=sptensor2spmat(J);

if isfield(model,'special_code')
    T=[T;model.moreT];
    J=[J;model.moreJ];
end


[mJ,nJ]=size(J);

if isempty(Wcols)
    Wcols=max(abs(J))';
    if min(Wcols)==0
        Wcols=ones(size(Wcols));
    end
    Wrows=max(abs([J*spdiags(1./Wcols,0,nJ,nJ)]'))';
    if min(Wrows)==0
        Wrows=ones(size(Wrows));
    end
end

scaledT=T./Wrows;
scaledJ=spdiags(1./Wrows,0,mJ,mJ)*J*spdiags(1./Wcols,0,nJ,nJ);



end

