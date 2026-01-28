function [ fxxxxxc,ind ] = chain5c_tensor(fv,fvv,fvvv,fvvvv,fvvvvv,vx,vxxc,vxxxc,vxxxxc,vxxxxxc...
    ,M1,M2,M3,M4,M5,M6,ind,n_ind,maxload,sum,varargin)
%chain5c is compressed.
%
% © Copyright, Oren Levintal, January 18, 2017.

if ~isempty(varargin)
    convertind=varargin{1};
    maxind=varargin{2};
    doi=1;
else
    doi=0;
end

if isempty(ind)
    ind=cell(14,1);
end

n_f=fv.tsize(1);
n_v=fv.tsize(2);
n_x=vx.tsize(2);
%fvvvvv*kron(vx,vx,vx,vx,vx)*U5

%fvvvvv*kron(kron(vx,vx,vx,vx)U4W4,vx))*U5
%fvvvvv*kron(kron(vx,vx,vx,vx)U4,vx))*kron(W4,I)*U5
fv_vvvv=col2ptr(fvvvvv,1); %n_vn_f,n_v,n_v,n_v,n_v
fv_vvvv=ptr1d(fv_vvvv); %n_vn_f,n_v,n_v,n_v,n_v
if doi==0
    [term1,ind{1}]=AkronBU4(fv_vvvv,vx,ind{1},n_ind,maxload,sum); %n_vn_f,unique4
else
    [term1,ind{1}]=AkronBU4i(fv_vvvv,vx,ind{1},n_ind,maxload,sum,reshape(repmat(convertind(:)',n_v,1),n_v*n_f,n_v),maxind); %n_vn_f,unique4
end
term1.ptr2d=fv_vvvv.ptr2d; 
term1=ptr2d(term1,n_f,n_v); 
term1=ptr2col(term1,2); %n_f,unique4,n_v
term1=col2ptr(term1,1); %unique4 n_f,n_v
term1=ptr1d(term1);
[term1b,ind{2}]=contraction1(term1,vx,ind{2},n_ind,maxload,sum); %unique4 n_f,n_x
term1b.ptr2d=term1.ptr2d;
term1b=ptr2d(term1b,term1b.ptr2d(1),term1b.ptr2d(2));
term1=ptr2col(term1b,2); % n_f,n_x,unique4
clear term1b
[term1,ind{3}]=contraction1(unfold(term1),M1,ind{3},n_ind,maxload,sum); %n_f,unique5

    % test
%     fvvvvv1=fvvvvv;
%     fvvvvv1.vals=fvvvvv1.vals(1,:);
%     vx1=vx;
%     vx1.vals=vx1.vals(1,:);
%     sfvvvvv=sptensor2spmat(unfold(fvvvvv1));
%     svx=sptensor2spmat(vx1);
%     sfv_vvvv=reshape(sfvvvvv,n_f*n_v,[]);
%     n_x=vx.tsize(2);
%     [U4,W4]=create_UW(n_x,4);
%     temp1=sfv_vvvv*kron(kron(svx,svx),kron(svx,svx))*U4;
%     temp1=reshape(full(temp1),n_f,n_v,[]);
%     temp1=permute(temp1,[2,1,3]);
%     temp1=sparse(reshape(temp1,n_v*n_f,[])); %n_v,n_f,unique4
%     
%     temp1=reshape(temp1,n_v,[])'*svx; % n_f,unique4,n_x
%     unique4=nchoosek(n_x+3,4);
%     temp1=reshape(full(temp1),n_f,unique4,n_x);
%     temp1=permute(temp1,[1,3,2]);
%     temp1=sparse(reshape(temp1,n_f,unique4*n_x));
% U5=create_UW(n_x,5);
%     term1_1=term1;
%     term1_1.vals=term1_1.vals(1,:);
%     temp1=sfvvvvv*kron(svx,kron(kron(svx,svx),kron(svx,svx)))*U5;
%     oren=sptensor2spmat(unfold(term1_1))-temp1;
%     max(abs(oren(:)))


%fvvvv*kron(kron(vx,vx,vx)*U3,vxxc)*kron(W3,W2)*OMEGA5*U5
[term2,ind{4}]=AkronB1U3B2(fvvvv,vx,vxxc,ind{4},n_ind,maxload,sum);
[term2,ind{5}]=contraction1(term2,M2,ind{5},n_ind,maxload,sum);

% test

% fvvvv1=fvvvv;
% fvvvv1.vals=fvvvv1.vals(1,:);
% vx1=vx;
% vx1.vals=vx1.vals(1,:);
% vxxc1=vxxc;
% vxxc1.vals=vxxc1.vals(1,:);
% sfvvvv=sptensor2spmat(unfold(fvvvv1));
% svx=sptensor2spmat(vx1);
% svxxc=sptensor2spmat(vxxc1);
% 
% n_x=vx.tsize(2);
% [U2,W2]=create_UW(n_x,2);
% [U5]=create_UW(n_x,5);
% 
% OMEGA=create_OMEGA(n_x,5);
% temp2=sfvvvv*kron(kron(svx,svx),kron(svx,svxxc*W2))*OMEGA.OMEGA5*U5;
% temp2=sparse(temp2);
% term2_1=term2;
% term2_1.vals=term2_1.vals(1,:);
% oren=sptensor2spmat(unfold(term2_1))-temp2;
% max(abs(oren(:)))

sumterm=tplus(term1,term2); clear term1 term2

%fvvv*kron(kron(vx,vx)*U2,vxxxc)*kron(W2,W3)*OMEGA6*U5
[term3,ind{6}]=AkronB1U2B2(fvvv,vx,vxxxc,ind{6},n_ind,maxload,sum);
[term3,ind{7}]=contraction1(term3,M3,ind{7},n_ind,maxload,sum);

% test

% fvvv1=fvvv;
% fvvv1.vals=fvvv1.vals(1,:);
% vx1=vx;
% vx1.vals=vx1.vals(1,:);
% vxxxc1=vxxxc;
% vxxxc1.vals=vxxxc1.vals(1,:);
% sfvvv=sptensor2spmat(unfold(fvvv1));
% svx=sptensor2spmat(vx1);
% svxxxc=sptensor2spmat(vxxxc1);
% 
% n_x=vx.tsize(2);
% [U2,W2]=create_UW(n_x,2);
% [U3,W3]=create_UW(n_x,3);
% [U5]=create_UW(n_x,5);
% 
% OMEGA=create_OMEGA(n_x,5);
% temp3=sfvvv*kron(kron(svx,svx),svxxxc*W3)*OMEGA.OMEGA6*U5;
% temp3=sparse(temp3);
% term3_1=term3;
% term3_1.vals=term3_1.vals(1,:);
% oren=sptensor2spmat(unfold(term3_1))-temp3;
% max(abs(oren(:)))

sumterm=tplus(sumterm,term3); clear term3

%fvvv*kron(vx,kron(vxxc*W2,vxxc*W2))*OMEGA7*U5
%fvvv*kron(vx,kron(vxxc,vxxc))*kron(I,W2,W2)*OMEGA7*U5
%fvvv*kron(vx,kron(vxxc,vxxc)UW)*kron(I,W2,W2)*OMEGA7*U5
%fvvv*kron(vx,kron(vxxc,vxxc)U)*kron(I,W)*kron(I,W2,W2)*OMEGA7*U5
[term4,ind{8}]=AkronB1U2B2(fvvv,vxxc,vx,ind{8},n_ind,maxload,sum);
term4=permutecols(fold(term4,n_x,term4.tsize(2)/n_x),[2,1]);
[term4,ind{9}]=contraction1(unfold(term4),M4,ind{9},n_ind,maxload,sum);

% test

% fvvv1=fvvv;
% fvvv1.vals=fvvv1.vals(1,:);
% vx1=vx;
% vx1.vals=vx1.vals(1,:);
% vxxc1=vxxc;
% vxxc1.vals=vxxc1.vals(1,:);
% sfvvv=sptensor2spmat(unfold(fvvv1));
% svx=sptensor2spmat(vx1);
% svxxc=sptensor2spmat(vxxc1);
% 
% n_x=vx.tsize(2);
% [U2,W2]=create_UW(n_x,2);
% [U3,W3]=create_UW(n_x,3);
% [U5]=create_UW(n_x,5);
% 
% OMEGA=create_OMEGA(n_x,5);
% temp4=sfvvv*kron(svx,kron(svxxc*W2,svxxc*W2))*OMEGA.OMEGA7*U5;
% temp4=sparse(temp4);
% term4_1=term4;
% term4_1.vals=term4_1.vals(1,:);
% oren=sptensor2spmat(unfold(term4_1))-temp4;
% max(abs(oren(:)))

sumterm=tplus(sumterm,term4); clear term4

%fvv*kron(vx,vxxxxc)*kron(I,W4)*OMEGA8*U5
[term5,ind{10}]=contraction2(fvv,vx,vxxxxc,ind{10},n_ind,maxload,sum);
[term5,ind{11}]=contraction1(unfold(term5),M5,ind{11},n_ind,maxload,sum);
sumterm=tplus(sumterm,term5); clear term5

%fvv*kron(vxxc,vxxxc)*kron(W2,W3)*OMEGA9*U5
[term6,ind{12}]=contraction2(fvv,vxxc,vxxxc,ind{12},n_ind,maxload,sum);
[term6,ind{13}]=contraction1(unfold(term6),M6,ind{13},n_ind,maxload,sum);
sumterm=tplus(sumterm,term6); clear term6

%fv*vxxxxxc
[term7,ind{14}]=contraction1(fv,vxxxxxc,ind{14},n_ind,maxload,sum);

fxxxxxc=tplus(sumterm,term7);

if isfield(fv,'ptr2d')
    fxxxxxc.ptr2d=fv.ptr2d;
end

end

