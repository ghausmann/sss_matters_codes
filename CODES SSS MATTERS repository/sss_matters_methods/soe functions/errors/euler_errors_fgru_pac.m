function y = euler_errors_fgru_pac(P2,nxss,nyss,derivs,x_1,eta,epsi_nodes,weight_nodes,approx)
%This function computes Euler equation errors of the SOE (PAC model)
%It assumes the model of interest. 
%Expectations are approximated with monomials
%Inputs are P2 (vector of parameters),
%nxss (DSS of states in the auxiliary model), %nyss (DSS of controls in the
%auxiliary model), derivs (matrices of derivatives of the policy
%functions), x_1 (current states), eta matrix, epsi_nodes (nodes of monomials),
%weight_nodes (probabilities of each node), and approx (order or
%approximation).
%The function returns unit-free Euler errors at the current state of the
%economy
betta = P2(1);
nu = P2(2);
phi_d = P2(3);
D_bar = P2(4);

st = x_1-nxss;

%Evaluate decision rules for states and controls
mht = dr_ht(derivs,nxss,approx,st);
mgt = dr_gt(derivs,nyss,approx,st);

rt = mgt(3); %interest rate
ct = exp(mgt(1)); %consumption 
Dt1 = mht(1); %debt

my_mht = mht + eta*epsi_nodes';

%Evaluate future controls for each node of the monomials
my_gt1 = dr_gt(derivs,nyss,approx,(my_mht-nxss));
ct1 = exp(my_gt1(1,:));

%Compute expectations of the Euler equation
myexp_int_a = ct1.^-nu;
Et = weight_nodes'*(myexp_int_a');
%Return unit-free Euler errors
myexp_a = ( (1/(1+rt) - phi_d*( Dt1 - D_bar) )^-1 )*betta*Et;
y = abs( (myexp_a^(-1/nu))*(ct^-1) -1  );
   
