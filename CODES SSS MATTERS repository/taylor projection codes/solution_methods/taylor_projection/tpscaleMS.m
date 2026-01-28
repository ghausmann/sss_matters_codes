function [ scaledT,scaledJ,Wcols,Wrows,model ] = tpscaleMS( weighted_coeffs,x0,model,params,msparams,transition,c0,nep,P,Wcols,Wrows)
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

[T,J,model]=tpMS(coeffs,x0,model,params,msparams,transition,c0,nep,P);

if isempty(Wcols)
    Wcols=max(abs(J))';
    if min(Wcols)==0
        Wcols=ones(size(Wcols));
    end
    Wrows=max(abs([J*spdiags(1./Wcols,0,model.n_theta*model.n_regimes,model.n_theta*model.n_regimes)]'))';
    if min(Wrows)==0
        Wrows=ones(size(Wrows));
    end
end

scaledT=T./Wrows;
[mJ,nJ]=size(J);
% scaledJ=spdiags(1./Wrows,0,model.n_theta*model.n_regimes,model.n_theta*model.n_regimes)*J*spdiags(1./Wcols,0,model.n_theta*model.n_regimes,model.n_theta*model.n_regimes);
scaledJ=spdiags(1./Wrows,0,mJ,mJ)*J*spdiags(1./Wcols,0,nJ,nJ);


end

