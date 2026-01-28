phiv_theta_vecT=takerows(v_theta,Phivars);
phiv_theta_vecT.vals=reshape(repmat(phiv_theta_vecT.vals(:)',n_nodes,1),n_nodes*n_s,[]);

[ Phitheta,model.hind{hindi} ] = chain0_theta_tensor( Phiv_vec,phiv_theta_vecT,model.hind{hindi},n_ind,maxload,'vec' );
hindi=hindi+1;

if order>=1
    phivxT=takerows(vx,Phivars);
    phivxT.vals=reshape(repmat(phivxT.vals(:)',n_nodes,1),n_nodes*n_s,[]);

    phivx_theta_vecT=takerows(vx_theta,Phivars);
    phivx_theta_vecT=unfold(phivx_theta_vecT);
    phivx_theta_vecT.vals=reshape(repmat(phivx_theta_vecT.vals(:)',n_nodes,1),n_nodes*n_s,[]);
    
    [ Phix,model.hind{hindi} ] = chain1_tensor( Phiv_vec,phivxT,model.hind{hindi},n_ind,maxload,'vec' );
    hindi=hindi+1;

    [ Phixtheta,model.hind{hindi} ] = chain1_theta_tensor( Phiv_vec,Phivv_vec,...
        phivxT,phiv_theta_vecT,phivx_theta_vecT,...
        model.hind{hindi},n_ind,maxload,'vec' );
    hindi=hindi+1;

end
if order>=2
    phivxxcT=takerows(vxxc,Phivars);
    phivxxcT.vals=reshape(repmat(phivxxcT.vals(:)',n_nodes,1),n_nodes*n_s,[]);
    
    phivxxc_theta_vecT=takerows(vxxc_theta,Phivars);
    phivxxc_theta_vecT=unfold(phivxxc_theta_vecT);
    phivxxc_theta_vecT.vals=reshape(repmat(phivxxc_theta_vecT.vals(:)',n_nodes,1),n_nodes*n_s,[]);

    [phivxT,model.hind4{hindi4}]=colsort(phivxT,model.hind4{hindi4});
    hindi4=hindi4+1;
    [ Phixxc,model.hind{hindi} ] = chain2c_tensor( Phiv_vec,Phivv_vec,...
        phivxT,phivxxcT,model.hind{hindi},n_ind,maxload,'vec' );
    hindi=hindi+1;

    [ Phixxctheta,model.hind{hindi} ] = chain2c_theta_tensor( Phiv_vec,Phivv_vec,Phivvv_vec,...
        phivxT,phivxxcT,...
        phiv_theta_vecT,phivx_theta_vecT,phivxxc_theta_vecT,...
        model.chain2c_theta_M2,...
        model.hind{hindi},n_ind,maxload,'vec' );
    hindi=hindi+1;

end
if order>=3
    phivxxxcT=takerows(vxxxc,Phivars);
    phivxxxcT.vals=reshape(repmat(phivxxxcT.vals(:)',n_nodes,1),n_nodes*n_s,[]);

    phivxxxc_theta_vecT=takerows(vxxxc_theta,Phivars);
    phivxxxc_theta_vecT=unfold(phivxxxc_theta_vecT);
    phivxxxc_theta_vecT.vals=reshape(repmat(phivxxxc_theta_vecT.vals(:)',n_nodes,1),n_nodes*n_s,[]);

    [ Phixxxc,model.hind{hindi} ] = chain3c_tensor( Phiv_vec,Phivv_vec,Phivvv_vec,...
        phivxT,phivxxcT,phivxxxcT,...
        model.x_chain3c_M2,...
        model.hind{hindi},n_ind,maxload,'vec' );
    hindi=hindi+1;

    [ Phixxxctheta,model.hind{hindi} ] = chain3c_theta_tensor( Phiv_vec,Phivv_vec,Phivvv_vec,Phivvvv_vec,...
        phivxT,phivxxcT,phivxxxcT,...
        phiv_theta_vecT,phivx_theta_vecT,phivxxc_theta_vecT,phivxxxc_theta_vecT,...
        model.chain3c_theta_M2,model.chain3c_theta_M3,model.chain3c_theta_M4,model.chain3c_theta_M5,...
        model.hind{hindi},n_ind,maxload,'vec' );
    hindi=hindi+1;

end
if order>=4
    phivxxxxcT=takerows(vxxxxc,Phivars);
    phivxxxxcT.vals=reshape(repmat(phivxxxxcT.vals(:)',n_nodes,1),n_nodes*n_s,[]);
 
    phivxxxxc_theta_vecT=takerows(vxxxxc_theta,Phivars);
    phivxxxxc_theta_vecT=unfold(phivxxxxc_theta_vecT);
    phivxxxxc_theta_vecT.vals=reshape(repmat(phivxxxxc_theta_vecT.vals(:)',n_nodes,1),n_nodes*n_s,[]);

    [phivxxcT,model.hind4{hindi4}]=colsort(phivxxcT,model.hind4{hindi4});
    hindi4=hindi4+1;        
    [ Phixxxxc,model.hind4{hindi4} ] = chain4c_tensor( Phiv_vec,Phivv_vec,Phivvv_vec,Phivvvv_vec,...
        phivxT,phivxxcT,phivxxxcT,phivxxxxcT,...
        model.x_chain4c_M2,model.x_chain4c_M3,model.x_chain4c_M4,...
        model.hind4{hindi4},n_ind,maxload,'vec' );
    hindi4=hindi4+1;

    
    [ Phixxxxctheta,model.hind4{hindi4} ] = chain4c_theta_tensor( Phiv_vec,Phivv_vec,Phivvv_vec,Phivvvv_vec,Phivvvvv_vec,...
        phivxT,phivxxcT,phivxxxcT,phivxxxxcT,...
        phiv_theta_vecT,phivx_theta_vecT,phivxxc_theta_vecT,phivxxxc_theta_vecT,phivxxxxc_theta_vecT,...
        model.chain4c_theta_M2,model.chain4c_theta_M3,model.chain4c_theta_M5,model.chain4c_theta_M6,model.chain4c_theta_M9,model.chain4c_theta_M10,...
        model.hind4{hindi4},n_ind,maxload,'vec' );
    hindi4=hindi4+1;

end
tothindi=tothindi+7;