function [T,J,model]=tpMS(coeffs,x0,model,params,msparams,transition,c0,nep,P)

if ~isequal(x0,c0)
    error('x0 must be equal to c0')
end

% msparams should have n_regimes columns

if size(msparams,2)~=model.n_regimes
    error('msparams should have n_regimes columns')
end

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

% ms conditions
if isfield(model,'jacobian')
    model.jacobian=model.jacobian;
end
[T,J,model]=tpvec(coeffsvec,x0vec,model,paramsvec,c0vec,nep,P,coeffspvec);


T.vals=reshape(T.vals,model.n_regimes,model.n_regimes,[]);
T.vals=T.vals.*reshape(repmat(transition,1,size(T.vals,3)),model.n_regimes,model.n_regimes,size(T.vals,3)); % multiply by the transition matrix

J.vals=reshape(J.vals,model.n_regimes,model.n_regimes,[]);
J.vals=J.vals.*reshape(repmat(transition,1,size(J.vals,3)),model.n_regimes,model.n_regimes,size(J.vals,3)); % multiply by the transition matrix
J.vals=reshape(J.vals,model.n_regimes^2,[]);

% sum across sp
T.vals=reshape(sum(T.vals,2),model.n_regimes,[]);

[rows,cols,vals]=tfind(T);
rows=double(rows);
cols=double(cols);

ncond=prod(T.tsize);
addrows=repmat((0:model.n_regimes-1)*ncond,length(rows),1);
newrows=repmat(rows,model.n_regimes,1)+addrows(:);
vals=vals';
vals=vals(:);

Trows=newrows(:);
Tcols=ones(length(Trows),1);
Tvals=vals;
Trowdim=model.n_regimes*ncond;

[rows,cols,vals]=tfind(J);

rows=double(rows);
cols=double(cols);

% Comment: I do not sum across sp. This will be done by the sparse
% function.

% find cols that belong to current period coeffs and future period coeffs
current_cols=logical(cols<=model.n_theta);
future_cols=logical(cols>model.n_theta);

newcols=repmat(cols(:)',model.n_regimes^2,1);

newcols(:,current_cols)=newcols(:,current_cols)+repmat((s(:)-1)*model.n_theta,1,nnz(current_cols));
newcols(:,future_cols)=newcols(:,future_cols)-model.n_theta+repmat((sp(:)-1)*model.n_theta,1,nnz(future_cols));

newrows=repmat(rows(:)',model.n_regimes^2,1);

newrows=reshape(newrows,model.n_regimes,[]);
newrows=newrows+repmat((0:model.n_regimes-1)'*ncond,1,size(newrows,2));

vals=vals(:);

Jrows=newrows(:);
Jcols=newcols(:);
Jvals=vals;

Jrowdim=model.n_regimes*ncond;
Jcoldim=model.n_regimes*model.n_theta;

T=sparse(Trows,Tcols,Tvals,Trowdim,1);
J=sparse(Jrows,Jcols,Jvals,Jrowdim,Jcoldim);

end