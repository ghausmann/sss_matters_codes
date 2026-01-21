%This function computes the DSS of capital in the original FGRU model
function Kstar = my_Kstar(my_dss_params)

alppha = my_dss_params(1);
r_bar = my_dss_params(2);
delta = my_dss_params(3);
etapar = my_dss_params(4);

%DSS of the model
Hstar = ( (1-alppha)*(alppha/((r_bar)+(delta)))^(alppha/(1-alppha))  )^(1/(etapar-1));
Kstar = Hstar*(alppha/((r_bar)+(delta)))^(1/(1-alppha));

