function [ M ] = chainsM( n,order,varargin )
%calculate coefficient matrices for compressed chain rules.
%
% © Copyright, Oren Levintal, June 13, 2016.

compress=1;
if ~isempty(varargin)
    compress=0;
end

if order==3
    M=cell(2,1);
    [~,W2]=create_UW(n,2);
    [U3,~]=create_UW(n,3);
    OMEGA=create_OMEGA(n,3);
    M{2}=spmat2sptensor(kron(speye(n),W2)*OMEGA.OMEGA1*U3);
elseif order==4
    M=cell(4,1);
    [~,W2]=create_UW(n,2);
    [~,W3]=create_UW(n,3);
    [U4,~]=create_UW(n,4);
    OMEGA=create_OMEGA(n,4);
    if compress==1
        M{2}=spmat2sptensor(kron(W2,W2)*(OMEGA.OMEGA2*U4)); OMEGA=rmfield(OMEGA,'OMEGA2');
        M{3}=spmat2sptensor(kron(speye(n),W3)*(OMEGA.OMEGA3*U4)); OMEGA=rmfield(OMEGA,'OMEGA3');
        [~,tempW]=create_UW(nchoosek(n+1,2),2);
        M{4}=spmat2sptensor(tempW*kron(W2,W2)*(OMEGA.OMEGA4*U4)); 
    elseif compress==0
        M{2}=spmat2sptensor(kron(W2,W2)*(OMEGA.OMEGA2)); OMEGA=rmfield(OMEGA,'OMEGA2');
        M{3}=spmat2sptensor(kron(speye(n),W3)*(OMEGA.OMEGA3)); OMEGA=rmfield(OMEGA,'OMEGA3');
        [~,tempW]=create_UW(nchoosek(n+1,2),2);
        M{4}=spmat2sptensor(tempW*kron(W2,W2)*(OMEGA.OMEGA4)); 
    end
elseif order==5
    M=cell(6,1);
    [~,W2]=create_UW(n,2);
    [~,W3]=create_UW(n,3);
    [U4,W4]=create_UW(n,4);
    [U5,W5]=create_UW(n,5);
    OMEGA=create_OMEGA(n,5);
    if compress==1
        M{1}=spmat2sptensor(kron(W4,speye(n))*U5);
        M{2}=spmat2sptensor(kron(W3,W2)*OMEGA.OMEGA5*U5); OMEGA=rmfield(OMEGA,'OMEGA5');
        M{3}=spmat2sptensor(kron(W2,W3)*OMEGA.OMEGA6*U5); OMEGA=rmfield(OMEGA,'OMEGA6');
        [~,tempW]=create_UW(nchoosek(n+1,2),2);
        M{4}=spmat2sptensor(kron(speye(n),tempW)*kron(kron(speye(n),W2),W2)*OMEGA.OMEGA7*U5); OMEGA=rmfield(OMEGA,'OMEGA7');
        M{5}=spmat2sptensor(kron(speye(n),W4)*OMEGA.OMEGA8*U5); OMEGA=rmfield(OMEGA,'OMEGA8');
        M{6}=spmat2sptensor(kron(W2,W3)*OMEGA.OMEGA9*U5); OMEGA=rmfield(OMEGA,'OMEGA9');

        %     elseif compress==0
%         M{2}=spmat2sptensor(kron(W2,W2)*(OMEGA.OMEGA2)); OMEGA=rmfield(OMEGA,'OMEGA2');
%         M{3}=spmat2sptensor(kron(speye(n),W3)*(OMEGA.OMEGA3)); OMEGA=rmfield(OMEGA,'OMEGA3');
%         [~,tempW]=create_UW(nchoosek(n+1,2),2);
%         M{4}=spmat2sptensor(tempW*kron(W2,W2)*(OMEGA.OMEGA4)); 
    end
else
    error('order is 3 or 4 or 5');
end

end

