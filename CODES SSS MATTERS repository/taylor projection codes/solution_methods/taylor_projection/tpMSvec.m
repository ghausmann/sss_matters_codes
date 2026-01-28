function [T,J,model]=tpMSvec(coeffs,x0,model,params,msparams,transition,c0,nep,P)

% This function is a vectorized version of tpMS. The MS system is computed
% for n_s different values of coeffs/x0/params, column-wise. The MS pattern
% is preserved (i.e. future coeffs are separate from current coeffs)
% INPUTS:
% coeffs (n_s,n_theta,n_regimes)
% x0 (n_s,n_x)
% params (n_s,n_params)
% c0 (n_s,n_x)
% msparams (n_msparams,n_regimes)
% OUTPUTS T and J are sparse tensors.
% (c) Oren Levintal, September 29, 2017.

if numel(x0)==model.n_x && numel(coeffs)==model.n_theta*model.n_regimes && numel(params)==model.n_params
    n_s=1;
else % x0,coeffs,params have a row dimension of 1 or n_s
    temp=[size(x0,1),size(coeffs,1),size(params,1),size(c0,1)];
    n_s=max(temp);
    if min(temp)~=n_s
        error('incompatible s dimension')
    end
end

if ~isequal(x0,c0)
    error('x0 must be equal to c0')
end

% msparams should have n_regimes columns

if size(msparams,2)~=model.n_regimes
    error('msparams should have n_regimes columns')
end

msparams=msparams'; %n_regimes,n_msparams

paramsvec=reshape(repmat(params(:)',model.n_regimes^2,1),model.n_regimes^2*n_s,[]); % n_regimes*n_regimes,n_s*n_params

[s,sp]=ndgrid([1:model.n_regimes],[1:model.n_regimes]);
s=s(:);
sp=sp(:);
paramsvec=[paramsvec,repmat([msparams(s,:),msparams(sp,:)],n_s,1)]; 

paramsvec=[paramsvec,zeros(model.n_regimes^2*n_s,model.n_logicalparams)];

coeffs=reshape(full(coeffs),n_s*model.n_theta,model.n_regimes);
coeffsvec=coeffs(:,s);
coeffsvec=reshape(coeffsvec,n_s,model.n_theta,model.n_regimes^2);
coeffsvec=permute(coeffsvec,[3,1,2]);
coeffsvec=reshape(coeffsvec,model.n_regimes^2*n_s,[]);

coeffspvec=coeffs(:,sp);
coeffspvec=reshape(coeffspvec,n_s,model.n_theta,model.n_regimes^2);
coeffspvec=permute(coeffspvec,[3,1,2]);
coeffspvec=reshape(coeffspvec,model.n_regimes^2*n_s,[]);


x0vec=repmat(x0(:)',model.n_regimes^2,1); %n_regimes^2,n_s*n_x
x0vec=reshape(x0vec,model.n_regimes^2*n_s,[]); %n_regimes^2*n_s,n_x
c0vec=repmat(c0(:)',model.n_regimes^2,1); %n_regimes^2,n_s*n_x
c0vec=reshape(c0vec,model.n_regimes^2*n_s,[]); %n_regimes^2*n_s,n_x

nep=full(nep);
P=full(P);

% ms conditions
if isfield(model,'jacobian')
    model.jacobian=model.jacobian;
end

% row dimensions should be n_regimes*n_regimes*n_s

[T,J,model]=tpvec(coeffsvec,x0vec,model,paramsvec,c0vec,nep,P,coeffspvec);


T.vals=reshape(T.vals,model.n_regimes,model.n_regimes,[]);
T.vals=T.vals.*reshape(repmat(transition,1,size(T.vals,3)),model.n_regimes,model.n_regimes,size(T.vals,3)); % multiply by the transition matrix

J.vals=reshape(J.vals,model.n_regimes,model.n_regimes,[]);
J.vals=J.vals.*reshape(repmat(transition,1,size(J.vals,3)),model.n_regimes,model.n_regimes,size(J.vals,3)); % multiply by the transition matrix
J.vals=reshape(J.vals,model.n_regimes^2*n_s,[]);

% sum across sp and reshape
T.vals=reshape(sum(T.vals,2),model.n_regimes*n_s,[]);

[rows,cols,vals]=tfind(T);
rows=double(rows);
cols=double(cols);

ncond=prod(T.tsize);
addrows=repmat((0:model.n_regimes-1)*ncond,length(rows),1);
newrows=repmat(rows,model.n_regimes,1)+addrows(:);
vals=reshape(vals,model.n_regimes,n_s,[]);
vals=permute(vals,[2,3,1]); %n_s,n_vals,n_regimes
vals=reshape(vals,n_s,[]); %n_s,n_vals*n_regimes


Trows=newrows(:);
Tcols=ones(length(Trows),1);
Tvals=vals;
Trowdim=model.n_regimes*ncond;




[rows,cols,vals]=tfind(J);

rows=double(rows);
cols=double(cols);

% Comment: I do not sum across sp. This is done later by post-multiplying
% the crs Jacobian matrix (with duplicate entries) by the identity matrix.

% find cols that belong to current period coeffs and future period coeffs
current_cols=logical(cols<=model.n_theta);
future_cols=logical(cols>model.n_theta);

newcols=repmat(cols(:)',model.n_regimes^2,1);

newcols(:,current_cols)=newcols(:,current_cols)+repmat((s(:)-1)*model.n_theta,1,nnz(current_cols));
newcols(:,future_cols)=newcols(:,future_cols)-model.n_theta+repmat((sp(:)-1)*model.n_theta,1,nnz(future_cols));

newrows=repmat(rows(:)',model.n_regimes^2,1);

newrows=reshape(newrows,model.n_regimes,[]);
newrows=newrows+repmat((0:model.n_regimes-1)'*ncond,1,size(newrows,2));

vals=reshape(vals,model.n_regimes^2,n_s,[]);
vals=permute(vals,[2,1,3]);
vals=reshape(vals,n_s,[]);

Jrows=newrows(:);
Jcols=newcols(:);
Jvals=vals;

Jrowdim=model.n_regimes*ncond;
Jcoldim=model.n_regimes*model.n_theta;

T=sptensor(Trows,Tcols,Tvals,Trowdim,1);
J=sptensor(Jrows,Jcols,Jvals,Jrowdim,Jcoldim);
J=contraction1(J,spteye(Jcoldim)); % this function sums across entries with same row and column (like the function sparse where there are duplicate entries)
% T=sparse(Trows,Tcols,Tvals,Trowdim,1);
% J=sparse(Jrows,Jcols,Jvals,Jrowdim,Jcoldim);

end