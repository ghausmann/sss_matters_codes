function [varout1,varout2]=contraction4(A,B4,B3,B2,B1,varargin)
%C=contraction4(A,B4,B3,B2,B1) calculates C=A*kron(B4,B3,B2,B1).
%[C,ind]=contraction4(A,B4,B3,B2,B1,ind,1,maxload) calculates C using one precomputed
%index ind. If ind=[] the function returns ind. The option maxload controls
%the maximum number of flops performed in parallel by the simd. This is
%relevant only if A or B1,B2,B3,B4 represent sets of s tensors.
%[C,ind]=contraction4(A,B4,B3,B2,B1,ind,2,maxload) calculates C using two precomputed
%indices ind.
%[C,ind]=contraction4(A,B4,B3,B2,B1,ind,1,maxload,'sum') sums across s if A
%or B1,B2 represent sets of s tensors. This option is useful for calculating
%expected value of C.
%[C,ind]=contraction4(A,B4,B3,B2,B1,ind,2,maxload,'sum') is similar but uses two
%precomputed indices.
%ind=contraction4(A,B4,B3,B2,B1,1) returns a structure ind with one field
%to speed up calculation.
%ind=contraction4(A,B4,B3,B2,B1,2) returns a structure ind with two fields. Speedup
%is larger compared to one field, but requires more memory.

if nargin<5 || nargin>9
    error('wrong number of input arguments')
end

if A.tsize(2)~=B1.tsize(1) || A.tsize(3)~=B2.tsize(1) || A.tsize(4)~=B3.tsize(1)
    error('incompatible dimensions')
end
dosum=0;
if nargin==9
    if strcmp(varargin{4},'sum')
        dosum=1;
    elseif ~strcmp(varargin{4},'vec')
        error(['ninth argument must be charater ' '''' 'sum' '''' ' or ' '''' 'vec' ''''])
    end
end
if nargin==5
    type=0; % index is not available
    maxload=intarray(1);
end
if nargin==6
    if varargin{1}==1
        type=1;
    elseif varargin{1}==2
        type=2;
    else
        error('number of precomputed indices should be 1 or 2')
    end
end
calculate_ind=0;
if nargin>=7
    if varargin{2}==1
        type=3;
    elseif varargin{2}==2
        type=4;
    else
        error('number of precomputed indices should be 1 or 2')
    end
    if isempty(varargin{1})
        calculate_ind=1;
    else
        ind=varargin{1};
        if ~isa(ind,'struct')
            error('index variable incorrectly defined')
        end
        if ~isfield(ind,'IC_1')
            error('index variable incorrectly defined')
        end
        if varargin{2}==2
            if ~isfield(ind,'JC_1')
                error('index variable incorrectly defined')
            end
        end
    end
end
if nargin==7
    maxload=intarray(1);
elseif nargin>=7
    maxload=intarray(varargin{3});
end

l=A.tsize(1);
m1=A.tsize(2);
m2=A.tsize(3);
m3=A.tsize(4);
n1=B1.tsize(2);
n2=B2.tsize(2);
n3=B3.tsize(2);
n4=B4.tsize(2);
IA=A.ptr;
JA=A.cols;
IB1=B1.ptr;
JB1=B1.cols;
IB2=B2.ptr;
JB2=B2.cols;
IB3=B3.ptr;
JB3=B3.cols;
IB4=B4.ptr;
JB4=B4.cols;
if type==0 % calculate A*kron(B4,B3,B2,B1) without precomputed index
    varout1=sptensor;
    [IC_1,IC_2,IC_3,IC_4] = contraction4IC_mex(l,m1,m2,m3,n1,n2,n3,n4,IA,JA,IB1,JB1,IB2,JB2,IB3,JB3,IB4,JB4);
    JCrows_1=IC_1(end)-1;
    NA=A.vals;
    NB1=B1.vals;
    NB2=B2.vals;
    NB3=B3.vals;
    NB4=B4.vals;
    if ~isreal(NA) || ~isreal(NB1) || ~isreal(NB2) || ~isreal(NB3) || ~isreal(NB4)
        error('complex numbers not supported')
    end
    [varout1.vals,JC_1] = contraction4vecI_mex(l,m1,m2,m3,n1,n2,n3,n4,IA,JA,NA,...
         IB1,JB1,NB1,IB2,JB2,NB2,IB3,JB3,NB3,IB4,JB4,NB4,...
         IC_1,IC_2,IC_3,IC_4,JCrows_1,...
         maxload,dosum); 
    varout1.ptr=IC_1;
    varout1.cols=JC_1;
    varout1.tsize=[l,n1,n2,n3,n4];
elseif type==1 % compute IC
    [varout1.IC_1,varout1.IC_2,varout1.IC_3,varout1.IC_4] = contraction4IC_mex(l,m1,m2,m3,n1,n2,n3,n4,IA,JA,IB1,JB1,IB2,JB2,IB3,JB3,IB4,JB4);
elseif type==2 % compute IC and JC
    [IC_1,IC_2,IC_3,IC_4] = contraction4IC_mex(l,m1,m2,m3,n1,n2,n3,n4,IA,JA,IB1,JB1,IB2,JB2,IB3,JB3,IB4,JB4);
    JCrows_1=IC_1(end)-1;
    JCrows_2=IC_2(end)-1;
    JCrows_3=IC_3(end)-1;
    JCrows_4=IC_4(end)-1;
    [ varout1.JC_1,varout1.JC_2,varout1.JC_3,varout1.JC_4 ] = contraction4JC_mex( l,m1,m2,m3,n1,n2,n3,n4,IA,JA,IB1,JB1,IB2,JB2,IB3,JB3,IB4,JB4,JCrows_1,JCrows_2,JCrows_3,JCrows_4,IC_2,IC_3,IC_4);
    varout1.IC_1=IC_1;
    varout1.IC_2=IC_2;
    varout1.IC_3=IC_3;
    varout1.IC_4=IC_4;
elseif type==3 % compute A*kron(B4,B3,B2,B1) with one precomputed index IC
    varout1=sptensor;
    if calculate_ind==1
        [IC_1,IC_2,IC_3,IC_4] = contraction4IC_mex(l,m1,m2,m3,n1,n2,n3,n4,IA,JA,IB1,JB1,IB2,JB2,IB3,JB3,IB4,JB4);
    else
        IC_1=ind.IC_1;
        IC_2=ind.IC_2;
        IC_3=ind.IC_3;
        IC_4=ind.IC_4;
    end
    JCrows_1=IC_1(end)-1;
    NA=A.vals;
    NB1=B1.vals;
    NB2=B2.vals;
    NB3=B3.vals;
    NB4=B4.vals;
    if ~isreal(NA) || ~isreal(NB1) || ~isreal(NB2) || ~isreal(NB3) || ~isreal(NB4)
        error('complex numbers not supported')
    end
    [varout1.vals,JC_1] = contraction4vecI_mex(l,m1,m2,m3,n1,n2,n3,n4,IA,JA,NA,...
         IB1,JB1,NB1,IB2,JB2,NB2,IB3,JB3,NB3,IB4,JB4,NB4,...
         IC_1,IC_2,IC_3,IC_4,JCrows_1,...
         maxload,dosum); 
    varout1.ptr=IC_1;
    varout1.cols=JC_1;
    varout1.tsize=[l,n1,n2,n3,n4];
    varout2.IC_1=IC_1;
    varout2.IC_2=IC_2;
    varout2.IC_3=IC_3;
    varout2.IC_4=IC_4;
elseif type==4 % compute A*kron(B4,B3,B2,B1) with two precomputed indices IC,JC
    varout1=sptensor;
    if calculate_ind==1
        [IC_1,IC_2,IC_3,IC_4] = contraction4IC_mex(l,m1,m2,m3,n1,n2,n3,n4,IA,JA,IB1,JB1,IB2,JB2,IB3,JB3,IB4,JB4);
        JCrows_1=IC_1(end)-1;
        JCrows_2=IC_2(end)-1;
        JCrows_3=IC_3(end)-1;
        JCrows_4=IC_4(end)-1;
        [ JC_1,JC_2,JC_3,JC_4 ] = contraction4JC_mex( l,m1,m2,m3,n1,n2,n3,n4,IA,JA,IB1,JB1,IB2,JB2,IB3,JB3,IB4,JB4,JCrows_1,JCrows_2,JCrows_3,JCrows_4,IC_2,IC_3,IC_4);
    else
        IC_1=ind.IC_1;
        IC_2=ind.IC_2;
        IC_3=ind.IC_3;
        IC_4=ind.IC_4;
        JC_1=ind.JC_1;
        JC_2=ind.JC_2;
        JC_3=ind.JC_3;
        JC_4=ind.JC_4;
    end
    NA=A.vals;
    NB1=B1.vals;
    NB2=B2.vals;
    NB3=B3.vals;
    NB4=B4.vals;
    if ~isreal(NA) || ~isreal(NB1) || ~isreal(NB2) || ~isreal(NB3) || ~isreal(NB4)
        error('complex numbers not supported')
    end
    varout1.vals = contraction4vec_mex(l,m1,m2,m3,n1,n2,n3,n4,IA,JA,NA,...
         IB1,JB1,NB1,IB2,JB2,NB2,IB3,JB3,NB3,IB4,JB4,NB4,...
         IC_1,JC_1,IC_2,JC_2,IC_3,JC_3,IC_4,JC_4,...
         maxload,dosum); 
    varout1.ptr=IC_1;
    varout1.cols=JC_1;
    varout1.tsize=[l,n1,n2,n3,n4];
    varout2.IC_1=IC_1;
    varout2.IC_2=IC_2;
    varout2.IC_3=IC_3;
    varout2.IC_4=IC_4;
    varout2.JC_1=JC_1;
    varout2.JC_2=JC_2;
    varout2.JC_3=JC_3;
    varout2.JC_4=JC_4;
end
if isfield(A,'ptr2d') % if the pointer of A is 2D varout is also 2D
    varout1.ptr2d=A.ptr2d;
end
end