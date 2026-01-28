function [ R ,ind] = chain4c_theta_tensor(fv,fvv,fvvv,fvvvv,fvvvvv,...
    vx,vxxc,vxxxc,vxxxxc,vtheta,vxtheta,vxxctheta,vxxxctheta,vxxxxctheta,...
    M2,M3,M5,M6,M9,M10,ind,n_ind,maxload,sum )
%
% © Copyright, Oren Levintal, January 23, 2017.

if isempty(ind)
    ind=cell(21,1);
end

n_f=fv.tsize(1);
n_v=fv.tsize(2);
n_x=vx.tsize(2);
n_theta=vtheta.tsize(2);
unique2=nchoosek(n_x+1,2);
unique3=nchoosek(n_x+2,3);
unique4=nchoosek(n_x+3,4);

n_s=max([size(fv.vals,1),size(vx.vals,1),size(vtheta.vals,1)]);

vxtheta=unfold(vxtheta);
vxxctheta=unfold(vxxctheta);
vxxxctheta=unfold(vxxxctheta);
vxxxxctheta=unfold(vxxxxctheta);

%fvvvvv*kron(vtheta,kron(vx,vx,vx,vx)*U4)

if ~isempty(fvvvvv.vals)
    fv_vvvv=col2ptr(fvvvvv,1); %n_vn_f,n_v,n_v,n_v,n_v
    fv_vvvv=ptr1d(fv_vvvv); %n_vn_f,n_v,n_v,n_v,n_v
    [term1,ind{1}]=AkronBU4(fv_vvvv,vx,ind{1},n_ind,maxload,sum); % n_vn_f,unique4
    term1.ptr2d=fv_vvvv.ptr2d;
    term1=ptr2d(term1,n_f,n_v); 
    term1=ptr2col(term1,2); %n_f,unique4,n_v
    term1=col2ptr(term1,1); %unique4 n_f,n_v
    term1=ptr1d(term1);
    [term1,ind{2}]=contraction1(term1,vtheta,ind{2},n_ind,maxload,sum); %unique4 n_f,n_theta
    term1=ptr2d(term1,n_f,unique4);
else
    term1=sptensor(n_f,[unique4,n_theta],n_s);
    term1=col2ptr(term1,1);
end

[term2,ind{3}]=AkronB1U3B2(fvvvv,vx,vxtheta,ind{3},n_ind,maxload,sum); %n_f,n_x,n_theta,unique3

term2=fold(term2,n_x,n_theta,unique3);

term2=ptr1d(unfold(col2ptr(term2,2))); %n_thetan_f,n_x,unique3


[term2,ind{4}]=contraction1(term2,M2,ind{4},n_ind,maxload,sum);
term2=ptr2d(term2,n_f,n_theta);% tsize=[n_f*n_theta,unique4]
term2=ptr2col(term2,2);% tsize=[n_f,unique4,n_theta]
term2=col2ptr(term2,1);% tsize=[n_f*unique4,n_theta]

sumterm=tplus(term1,term2,maxload);


clear term1 term2

% term3

fv_vvv=col2ptr(fvvvv,1); %n_f*n_v,n_v,n_v,n_v
fv_vvv=ptr1d(fv_vvv); %n_f*n_v,n_v,n_v,n_v
[term3,ind{5}]=AkronB1U2B2(fv_vvv,vx,vxxc,ind{5},n_ind,maxload,sum); %n_f*n_v,unique2,unique2
term3.ptr2d=fv_vvv.ptr2d;
term3=ptr2d(term3,n_f,n_v); %n_f*n_v,unique2*unique2
term3=ptr2col(term3,1); %n_f,n_v,unique2*unique2

term3=fold(term3,n_v,unique2*unique2);
term3=col2ptr(term3,2);
term3=ptr1d(term3);
[term3,ind{6}]=contraction1(term3,vtheta,ind{6},n_ind,maxload,sum);
term3=fold(ptr2col(ptr2d(term3,n_f,unique2*unique2),2),n_theta,unique2,unique2);

term3=col2ptr(term3,1); %n_f*n_theta,unique2,unique2
term3=ptr1d(unfold(term3)); %n_f*n_theta,unique2*unique2
[term3,ind{7}]=contraction1(term3,M3,ind{7},n_ind,maxload,sum);
term3=ptr2d(term3,n_f,n_theta); %n_f*n_theta,unique4
term3=ptr2col(term3,2);%n_f,unique4,n_theta
term3=col2ptr(term3,1);%n_f*unique4,n_theta

sumterm=tplus(sumterm,term3,maxload);

clear term3

% term4 

[term4,ind{8}]=AkronB1U2B2(fvvv,vx,vxxctheta,ind{8},n_ind,maxload,sum); %n_f,unique2,n_theta,unique2
term4=fold(term4,unique2,n_theta,unique2);
term4=col2ptr(term4,2); %n_thetan_f,unique2,unique2
term4=ptr1d(unfold(term4)); %n_thetan_f,unique2*unique2
[term4,ind{9}]=contraction1(term4,M3,ind{9},n_ind,maxload,sum);
term4=ptr2d(term4,n_f,n_theta); %n_f*n_theta,unique4
term4=ptr2col(term4,2);%n_f,unique4,n_theta
term4=col2ptr(term4,1);%n_f*unique4,n_theta

sumterm=tplus(sumterm,term4,maxload);



clear term4

% term5

[term5,ind{10}]=contraction3(fvvv,vx,vxtheta,vxxc,ind{10},n_ind,maxload,sum); %n_f,unique2,n_x,n_theta,n_x
term5=fold(term5,unique2,n_x,n_theta,n_x);
term5=col2ptr(term5,3); %n_thetan_f,unique2,n_x,n_x
term5=ptr1d(unfold(term5)); %n_thetan_f,unique2*n_x*n_x
[term5,ind{11}]=contraction1(term5,M5,ind{11},n_ind,maxload,sum);
term5=ptr2d(term5,n_f,n_theta); %n_f*n_theta,unique4
term5=ptr2col(term5,2);%n_f,unique4,n_theta
term5=col2ptr(term5,1);%n_f*unique4,n_theta

sumterm=tplus(sumterm,term5,maxload);


clear term5

% term6

[term6,ind{12}]=contraction3(fvvv,vx,vxxxc,vtheta,ind{12},n_ind,maxload,sum); %n_f,n_theta,unique3,n_x
term6=fold(term6,n_theta,unique3,n_x);
term6=col2ptr(term6,1); %n_thetan_f,unique3,n_x


% term7

[term7,ind{13}]=contraction2(fvv,vx,vxxxctheta,ind{13},n_ind,maxload,sum); %n_f,unique3,n_theta,n_x
term7=fold(term7,unique3,n_theta,n_x);
term7=col2ptr(term7,2); %n_thetan_f,unique3,n_x

term6_7_8=tplus(unfold(term6),unfold(term7),maxload);

clear term6 term7

[term8,ind{14}]=contraction2(fvv,vxtheta,vxxxc,ind{14},n_ind,maxload,sum); %n_f,unique3,n_x,n_theta
term8=fold(term8,unique3,n_x,n_theta);
term8=col2ptr(term8,3); %n_thetan_f,unique3,n_x

term6_7_8=tplus(term6_7_8,unfold(term8),maxload); %n_thetan_f,unique3,n_x

clear term8

term6_7_8=unfold(term6_7_8); %n_thetan_f,unique3,n_x
term6_7_8=ptr1d(term6_7_8);
[term6_7_8,ind{15}]=contraction1(term6_7_8,M6,ind{15},n_ind,maxload,sum);
term6_7_8=ptr2d(term6_7_8,n_f,n_theta); %n_f*n_theta,unique4
term6_7_8=ptr2col(term6_7_8,2);%n_f,unique4,n_theta
term6_7_8=col2ptr(term6_7_8,1);%n_f*unique4,n_theta

sumterm=tplus(sumterm,term6_7_8,maxload);


clear term6_7_8

% term9
unique2new=nchoosek(unique2+1,2);
[term9,ind{16}]=AkronB1U2B2(fvvv,vxxc,vtheta,ind{16},n_ind,maxload,sum); %n_f,n_theta,unique2new
term9=fold(term9,n_theta,unique2new);
term9=col2ptr(term9,1); %n_thetan_f,unique2new
term9=ptr1d(unfold(term9)); %n_thetan_f,unique2new
[term9,ind{17}]=contraction1(term9,M9,ind{17},n_ind,maxload,sum);
term9=ptr2d(term9,n_f,n_theta); %n_f*n_theta,unique4
term9=ptr2col(term9,2);%n_f,unique4,n_theta
term9=col2ptr(term9,1);%n_f*unique4,n_theta

sumterm=tplus(sumterm,term9,maxload);

clear term9

% term10

[term10,ind{18}]=contraction2(fvv,vxxc,vxxctheta,ind{18},n_ind,maxload,sum); %n_f,unique2,n_theta,unique2
term10=fold(term10,unique2,n_theta,unique2);
term10=col2ptr(term10,2); %n_thetan_f,unique2,unique2
term10=ptr1d(unfold(term10)); %n_thetan_f,unique2*unique2
[term10,ind{19}]=contraction1(term10,M10,ind{19},n_ind,maxload,sum);
term10=ptr2d(term10,n_f,n_theta); %n_f*n_theta,unique4
term10=ptr2col(term10,2);%n_f,unique4,n_theta
term10=col2ptr(term10,1);%n_f*unique4,n_theta

sumterm=tplus(sumterm,term10,maxload);


clear term10

% term11

[term11,ind{20}]=contraction2(fvv,vtheta,vxxxxc,ind{20},n_ind,maxload,sum); % n_f,unique4,n_theta
term11=col2ptr(term11,1);%n_f*unique4,n_theta

sumterm=tplus(sumterm,term11,maxload);


clear term11

% term12

[term12,ind{21}]=contraction1(fv,vxxxxctheta,ind{21},n_ind,maxload,sum); % n_f,unique4,n_theta
term12=fold(term12,unique4,n_theta);
term12=col2ptr(term12,1);%n_f*unique4,n_theta

R=tplus(sumterm,term12,maxload);

R=ptr2col(R,1);

end

