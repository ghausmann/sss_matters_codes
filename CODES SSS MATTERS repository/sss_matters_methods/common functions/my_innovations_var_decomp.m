function [y1,y2] = my_innovations_var_decomp(T,replications,P1,P2,vol_out,level_out,prod_out)
%This function generates the pseudo-random iid innovations for the variance decomposition.
%Note that (i) shocks affecting the interest rate are restricted to the [-1,1] interval; 
%and (ii) we only admits sequence that produce non-negative interest rates.

rho_eps_tb=P1(1);
sigma_tb_bar=P1(2);
rho_sigma_tb=P1(3);
eta_tb=P1(4);

rho_eps_r=P2(1);
sigma_r_bar=P2(2);
rho_sigma_r=P2(3);
eta_r=P2(4);
r_bar=P2(5);
k = P2(6);

Sigma = [1 0 k 0 0;0 1 0 0 0;k 0 1 0 0;0 0 0 1 0;0 0 0 0 1];
    
    
innovations = zeros(5,T,replications);
rates = zeros(replications,T);

s=1;
while s<replications+1
     shocks = randn(5,T);
%     keyboard;
    %shocks = mvnrnd(zeros(5,1),Sigma,T)';
%     shocks(shocks<-1)=-1;
%     shocks(shocks>1)=1;
    

    shocks_rates = shocks(1:4,:);
    shocks_rates(shocks_rates<-1)=-1;
    shocks_rates(shocks_rates>1)=1;
    shocks(1:4,:) = shocks_rates;
   
%     shocks(shocks(1:4,:)<-1)=-1;
%     shocks(shocks(1:4,:)>1)=1;
     %keyboard;
    u_sigma_r = shocks(1,:);
    u_sigma_tb = shocks(2,:);
    u_r = shocks(3,:);
    u_tb = shocks(4,:);
    
    
    sigma_r = zeros(1,T+1);
    sigma_tb = zeros(1,T+1);
    sigma_tb(1) = sigma_tb_bar;
    sigma_r(1) = sigma_r_bar;
    eps_r = zeros(1,T+1);
    eps_tb = zeros(1,T+1);
    r = r_bar + zeros(1,T+1);
    

    for t = 1:T

        sigma_r(t+1)=(1-rho_sigma_r)*sigma_r_bar+rho_sigma_r*sigma_r(t)+eta_r*u_sigma_r(t);
        sigma_tb(t+1)=(1-rho_sigma_tb)*sigma_tb_bar+rho_sigma_tb*sigma_tb(t)+eta_tb*u_sigma_tb(t);
        eps_r(t+1)=rho_eps_r*eps_r(t)+exp(sigma_r(t+1))*u_r(t);
        eps_tb(t+1)=rho_eps_tb*eps_tb(t)+exp(sigma_tb(t+1))*u_tb(t);
        r(t+1) = r(t+1) + eps_tb(t+1) + eps_r(t+1);

    end
    
    my_min_r = min(r);
    if my_min_r>0
       
        if vol_out==1
            shocks(1:2,:)=0;
        end
        if level_out==1
            shocks(3:4,:)=0;
        end
        if prod_out==1
            shocks(5,:)=0;
        end
        
        
        
        innovations(:,:,s) = shocks;
        rates(s,:) = r(1:end-1);
        s = s+1;
                
    end
    
end

y1 = innovations;
y2 = rates;
