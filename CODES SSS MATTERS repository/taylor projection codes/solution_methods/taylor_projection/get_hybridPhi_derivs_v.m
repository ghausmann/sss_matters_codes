% compute derivatives of Phi w.r.t v by high order chain rules.
n_Phivars=length(model.Phivars);
Phivars=model.Phivars;

count_Phi_n=model.Phi_n-1;
if count_Phi_n<0
    count_Phi_n=0;
end

if ~isfield(model,'hind')
    model.hind=cell(100+1+8+8+count_Phi_n*(1+3+4+6)+1+4+6+9,1);
    model.hind4=cell(100+count_Phi_n,1);
end


hindi=1;
hindi4=1;
tothindi=1;

% extract PhiPIv,PhiPivv,... from PIz,PIzz,...
if order>=0
    [PhiPIz_vec,model.hind{hindi}]=extract(PIz_vec,model.Phizvars,model.Phizvars,0,model.hind{hindi});
    hindi=hindi+1;
    [PhiPIv_vec,model.hind{hindi}]=extract(PIz_vec,model.Phizvars,model.Phivars,0,model.hind{hindi});
    hindi=hindi+1;  
    fname=['Phi' model.fname '_tilf_d1'];
    Phiz_vec=feval(fname,[nv_vec(:,Phivars),nu_vec(:,model.Phiuvars)],params_vec,model.Phi_ind_u);
end
if order>=1
    [PhiPIzz_vec,model.hind{hindi}]=extract(PIzz_vec,model.Phizvars,model.Phizvars,0,model.hind{hindi});
    hindi=hindi+1;
    [PhiPIvvc_vec,model.hind{hindi}]=extract(PIzz_vec,model.Phizvars,model.Phivars,1,model.hind{hindi});
    hindi=hindi+1;
    fname=['Phi' model.fname '_tilf_d2'];
    Phizz_vec=feval(fname,[nv_vec(:,Phivars),nu_vec(:,model.Phiuvars)],params_vec,model.Phi_ind_u);
end
if order>=2
    [PhiPIzzz_vec,model.hind{hindi}]=extract(PIzzz_vec,model.Phizvars,model.Phizvars,0,model.hind{hindi});
    hindi=hindi+1;
    [PhiPIvvvc_vec,model.hind{hindi}]=extract(PIzzz_vec,model.Phizvars,model.Phivars,1,model.hind{hindi});
    hindi=hindi+1;
    fname=['Phi' model.fname '_tilf_d3'];
    Phizzz_vec=feval(fname,[nv_vec(:,Phivars),nu_vec(:,model.Phiuvars)],params_vec,model.Phi_ind_u);
    
    if ~isfield(model,'Phivars_chain3c_M2')
        [ tempM ] = chainsM( n_Phivars,3 );
        model.Phivars_chain3c_M2=tempM{2};
        clear tempM
    end
end
if order>=3
    [PhiPIzzzz_vec,model.hind{hindi}]=extract(PIzzzz_vec,model.Phizvars,model.Phizvars,0,model.hind{hindi});
    hindi=hindi+1;
    [PhiPIvvvvc_vec,model.hind{hindi}]=extract(PIzzzz_vec,model.Phizvars,model.Phivars,1,model.hind{hindi});
    hindi=hindi+1;
    fname=['Phi' model.fname '_tilf_d4'];
    Phizzzz_vec=feval(fname,[nv_vec(:,Phivars),nu_vec(:,model.Phiuvars)],params_vec,model.Phi_ind_u);
    if ~isfield(model,'Phivars_chain4c_M2')
        [ tempM ] = chainsM( n_Phivars,4 );
        model.Phivars_chain4c_M2=tempM{2};
        model.Phivars_chain4c_M3=tempM{3};
        model.Phivars_chain4c_M4=tempM{4};
        clear tempM
    end
end
if order>=4
    [PhiPIzzzzz_vec,model.hind4{hindi4}]=extract(PIzzzzz_vec,model.Phizvars,model.Phizvars,0,model.hind4{hindi4});
    hindi4=hindi4+1;
    [PhiPIvvvvvc_vec,model.hind4{hindi4}]=extract(PIzzzzz_vec,model.Phizvars,model.Phivars,1,model.hind4{hindi4});
    hindi4=hindi4+1;
    fname=['Phi' model.fname '_tilf_d5'];
    Phizzzzz_vec=feval(fname,[nv_vec(:,Phivars),nu_vec(:,model.Phiuvars)],params_vec,model.Phi_ind_u);
end
tothindi=tothindi+8;
% hindi=tothindi;
    


if order==0
    hindi=tothindi;
    for i=2:model.Phi_n
        [PhiPIv_vec,model.hind{hindi}]=chain1_tensor(PhiPIz_vec,PhiPIv_vec,model.hind{hindi},n_ind,maxload,'vec');
        hindi=hindi+1;
    end
    [Phiv_vec,model.hind{hindi}]=chain1_tensor(Phiz_vec,PhiPIv_vec,model.hind{hindi},n_ind,maxload,'vec');
    hindi=hindi+1;

%     stochfv_vec.vals=multcol(stochfv_vec.vals,Pvec);
% 
%     model.stochfv_vec=stochfv_vec;

elseif order==1
    hindi=tothindi+count_Phi_n*(1)+1;
    for i=2:model.Phi_n
        [PhiPIv_vec,model.hind{hindi}]=colsort(PhiPIv_vec,model.hind{hindi});
        hindi=hindi+1;
        [PhiPIvvc_vec,model.hind{hindi}]=chain2c_tensor(PhiPIz_vec,PhiPIzz_vec,PhiPIv_vec,PhiPIvvc_vec,model.hind{hindi},n_ind,maxload,'vec');
        hindi=hindi+1;
        [PhiPIv_vec,model.hind{hindi}]=chain1_tensor(PhiPIz_vec,PhiPIv_vec,model.hind{hindi},n_ind,maxload,'vec');
        hindi=hindi+1;
    end
    [PhiPIv_vec,model.hind{hindi}]=colsort(PhiPIv_vec,model.hind{hindi});
    hindi=hindi+1;

    [Phiv_vec,model.hind{hindi}]=chain1_tensor(Phiz_vec,PhiPIv_vec,model.hind{hindi},n_ind,maxload,'vec');
    hindi=hindi+1;
    [Phivvc_vec,model.hind{hindi}]=chain2c_tensor(Phiz_vec,Phizz_vec,PhiPIv_vec,PhiPIvvc_vec,model.hind{hindi},n_ind,maxload,'vec');
    hindi=hindi+1;

%     stochfv_vec.vals=multcol(stochfv_vec.vals,Pvec);
%     stochfvvc_vec.vals=multcol(stochfvvc_vec.vals,Pvec);

    [ Phivv_vec,model.hind{hindi} ] = uncompressderivs( Phivvc_vec,2,n_Phivars,model.Phiv(:,Phivars),model.hind{hindi} );
    hindi=hindi+1;
    clear Phivvc_vec

elseif order==2
    hindi=tothindi+count_Phi_n*(1+3)+1+4;
    for i=2:model.Phi_n
        [PhiPIv_vec,model.hind{hindi}]=colsort(PhiPIv_vec,model.hind{hindi});
        hindi=hindi+1;
        [PhiPIvvvc_vec,model.hind{hindi}]=chain3c_tensor(PhiPIz_vec,PhiPIzz_vec,PhiPIzzz_vec,...
            PhiPIv_vec,PhiPIvvc_vec,PhiPIvvvc_vec,model.Phivars_chain3c_M2,model.hind{hindi},n_ind,maxload,'vec');
        hindi=hindi+1;
        [PhiPIvvc_vec,model.hind{hindi}]=chain2c_tensor(PhiPIz_vec,PhiPIzz_vec,PhiPIv_vec,PhiPIvvc_vec,model.hind{hindi},n_ind,maxload,'vec');
        hindi=hindi+1;
        [PhiPIv_vec,model.hind{hindi}]=chain1_tensor(PhiPIz_vec,PhiPIv_vec,model.hind{hindi},n_ind,maxload,'vec');
        hindi=hindi+1;
    end
    [PhiPIv_vec,model.hind{hindi}]=colsort(PhiPIv_vec,model.hind{hindi});
    hindi=hindi+1;

    [Phiv_vec,model.hind{hindi}]=chain1_tensor(Phiz_vec,PhiPIv_vec,model.hind{hindi},n_ind,maxload,'vec');
    hindi=hindi+1;
    [Phivvc_vec,model.hind{hindi}]=chain2c_tensor(Phiz_vec,Phizz_vec,PhiPIv_vec,PhiPIvvc_vec,model.hind{hindi},n_ind,maxload,'vec');
    hindi=hindi+1;
    [Phivvvc_vec,model.hind{hindi}]=chain3c_tensor(Phiz_vec,Phizz_vec,Phizzz_vec,...
        PhiPIv_vec,PhiPIvvc_vec,PhiPIvvvc_vec,...
        model.Phivars_chain3c_M2,...
        model.hind{hindi},n_ind,maxload,'vec');
    hindi=hindi+1;

%     stochfv_vec.vals=multcol(stochfv_vec.vals,Pvec);
%     stochfvvc_vec.vals=multcol(stochfvvc_vec.vals,Pvec);
%     stochfvvvc_vec.vals=multcol(stochfvvvc_vec.vals,Pvec);

    [ Phivv_vec,model.hind{hindi} ] = uncompressderivs( Phivvc_vec,2,n_Phivars,model.Phiv(:,Phivars),model.hind{hindi} );
    hindi=hindi+1;
    clear Phivvc_vec

    [ Phivvv_vec,model.hind{hindi} ] = uncompressderivs( Phivvvc_vec,3,n_Phivars,model.Phiv(:,Phivars),model.hind{hindi} );
    hindi=hindi+1;
    clear Phivvvc_vec

elseif order==3
    hindi=tothindi+count_Phi_n*(1+3+4)+1+4+6;
    for i=2:model.Phi_n
        [PhiPIv_vec,model.hind{hindi}]=colsort(PhiPIv_vec,model.hind{hindi});
        hindi=hindi+1;
        [PhiPIvvc_vec,model.hind{hindi}]=colsort(PhiPIvvc_vec,model.hind{hindi});
        hindi=hindi+1;
        [PhiPIvvvvc_vec,model.hind{hindi}]=chain4c_tensor(PhiPIz_vec,PhiPIzz_vec,PhiPIzzz_vec,PhiPIzzzz_vec,...
            PhiPIv_vec,PhiPIvvc_vec,PhiPIvvvc_vec,PhiPIvvvvc_vec,...
            model.Phivars_chain4c_M2,model.Phivars_chain4c_M3,model.Phivars_chain4c_M4,...
            model.hind{hindi},n_ind,maxload,'vec',model.Phizz,model.maxzz);
        hindi=hindi+1;
        [PhiPIvvvc_vec,model.hind{hindi}]=chain3c_tensor(PhiPIz_vec,PhiPIzz_vec,PhiPIzzz_vec,...
            PhiPIv_vec,PhiPIvvc_vec,PhiPIvvvc_vec,model.Phivars_chain3c_M2,model.hind{hindi},n_ind,maxload,'vec');
        hindi=hindi+1;
        [PhiPIvvc_vec,model.hind{hindi}]=chain2c_tensor(PhiPIz_vec,PhiPIzz_vec,PhiPIv_vec,PhiPIvvc_vec,model.hind{hindi},n_ind,maxload,'vec');
        hindi=hindi+1;
        [PhiPIv_vec,model.hind{hindi}]=chain1_tensor(PhiPIz_vec,PhiPIv_vec,model.hind{hindi},n_ind,maxload,'vec');
        hindi=hindi+1;
    end
    [PhiPIv_vec,model.hind{hindi}]=colsort(PhiPIv_vec,model.hind{hindi});
    hindi=hindi+1;
    [PhiPIvvc_vec,model.hind{hindi}]=colsort(PhiPIvvc_vec,model.hind{hindi});
    hindi=hindi+1;

    [Phiv_vec,model.hind{hindi}]=chain1_tensor(Phiz_vec,PhiPIv_vec,model.hind{hindi},n_ind,maxload,'vec');
    hindi=hindi+1;
    [Phivvc_vec,model.hind{hindi}]=chain2c_tensor(Phiz_vec,Phizz_vec,PhiPIv_vec,PhiPIvvc_vec,model.hind{hindi},n_ind,maxload,'vec');
    hindi=hindi+1;
    [Phivvvc_vec,model.hind{hindi}]=chain3c_tensor(Phiz_vec,Phizz_vec,Phizzz_vec,...
        PhiPIv_vec,PhiPIvvc_vec,PhiPIvvvc_vec,...
        model.Phivars_chain3c_M2,...
        model.hind{hindi},n_ind,maxload,'vec');
    hindi=hindi+1;

    [Phivvvvc_vec,model.hind{hindi}]=chain4c_tensor(Phiz_vec,Phizz_vec,Phizzz_vec,Phizzzz_vec,...
        PhiPIv_vec,PhiPIvvc_vec,PhiPIvvvc_vec,PhiPIvvvvc_vec,...
        model.Phivars_chain4c_M2,model.Phivars_chain4c_M3,model.Phivars_chain4c_M4,...
        model.hind{hindi},n_ind,maxload,'vec',model.tilPhiz,model.maxtilPhiz);
    hindi=hindi+1;

%     stochfv_vec.vals=multcol(stochfv_vec.vals,Pvec);
%     stochfvvc_vec.vals=multcol(stochfvvc_vec.vals,Pvec);
%     stochfvvvc_vec.vals=multcol(stochfvvvc_vec.vals,Pvec);
%     stochfvvvvc_vec.vals=multcol(stochfvvvvc_vec.vals,Pvec);

    [ Phivv_vec,model.hind{hindi} ] = uncompressderivs( Phivvc_vec,2,n_Phivars,model.Phiv(:,Phivars),model.hind{hindi} );
    hindi=hindi+1;
    clear Phivvc_vec

    [ Phivvv_vec,model.hind{hindi} ] = uncompressderivs( Phivvvc_vec,3,n_Phivars,model.Phiv(:,Phivars),model.hind{hindi} );
    hindi=hindi+1;
    clear Phivvvc_vec

    [ Phivvvv_vec,model.hind{hindi} ] = uncompressderivs( Phivvvvc_vec,4,n_Phivars,model.Phiv(:,Phivars),model.hind{hindi} );
    hindi=hindi+1;
    clear Phivvvvc_vec
    
elseif order==4
    hindi=tothindi+count_Phi_n*(1+3+4)+1+4+6;
    for i=2:model.Phi_n
        [PhiPIv_vec,model.hind{hindi}]=colsort(PhiPIv_vec,model.hind{hindi});
        hindi=hindi+1;
        [PhiPIvvc_vec,model.hind{hindi}]=colsort(PhiPIvvc_vec,model.hind{hindi});
        hindi=hindi+1;
        [PhiPIvvvvvc_vec,model.hind4{hindi4}]=chain5c_tensor(PhiPIz_vec,PhiPIzz_vec,PhiPIzzz_vec,PhiPIzzzz_vec,PhiPIzzzzz_vec,...
            PhiPIv_vec,PhiPIvvc_vec,PhiPIvvvc_vec,PhiPIvvvvc_vec,PhiPIvvvvvc_vec,...
            model.Phivars_chain5c_M1,model.Phivars_chain5c_M2,model.Phivars_chain5c_M3,model.Phivars_chain5c_M4,model.Phivars_chain5c_M5,model.Phivars_chain5c_M6,...
            model.hind4{hindi4},n_ind,maxload,'vec',model.Phizz,model.maxzz);
        hindi4=hindi4+1;
        [PhiPIvvvvc_vec,model.hind{hindi}]=chain4c_tensor(PhiPIz_vec,PhiPIzz_vec,PhiPIzzz_vec,PhiPIzzzz_vec,...
            PhiPIv_vec,PhiPIvvc_vec,PhiPIvvvc_vec,PhiPIvvvvc_vec,...
            model.Phivars_chain4c_M2,model.Phivars_chain4c_M3,model.Phivars_chain4c_M4,...
            model.hind{hindi},n_ind,maxload,'vec',model.Phizz,model.maxzz);
        hindi=hindi+1;
        [PhiPIvvvc_vec,model.hind{hindi}]=chain3c_tensor(PhiPIz_vec,PhiPIzz_vec,PhiPIzzz_vec,...
            PhiPIv_vec,PhiPIvvc_vec,PhiPIvvvc_vec,model.Phivars_chain3c_M2,model.hind{hindi},n_ind,maxload,'vec');
        hindi=hindi+1;
        [PhiPIvvc_vec,model.hind{hindi}]=chain2c_tensor(PhiPIz_vec,PhiPIzz_vec,PhiPIv_vec,PhiPIvvc_vec,model.hind{hindi},n_ind,maxload,'vec');
        hindi=hindi+1;
        [PhiPIv_vec,model.hind{hindi}]=chain1_tensor(PhiPIz_vec,PhiPIv_vec,model.hind{hindi},n_ind,maxload,'vec');
        hindi=hindi+1;
    end
    [PhiPIv_vec,model.hind{hindi}]=colsort(PhiPIv_vec,model.hind{hindi});
    hindi=hindi+1;
    [PhiPIvvc_vec,model.hind{hindi}]=colsort(PhiPIvvc_vec,model.hind{hindi});
    hindi=hindi+1;

    [Phiv_vec,model.hind{hindi}]=chain1_tensor(Phiz_vec,PhiPIv_vec,model.hind{hindi},n_ind,maxload,'vec');
    hindi=hindi+1;
    [Phivvc_vec,model.hind{hindi}]=chain2c_tensor(Phiz_vec,Phizz_vec,PhiPIv_vec,PhiPIvvc_vec,model.hind{hindi},n_ind,maxload,'vec');
    hindi=hindi+1;
    [Phivvvc_vec,model.hind{hindi}]=chain3c_tensor(Phiz_vec,Phizz_vec,Phizzz_vec,...
        PhiPIv_vec,PhiPIvvc_vec,PhiPIvvvc_vec,...
        model.Phivars_chain3c_M2,...
        model.hind{hindi},n_ind,maxload,'vec');
    hindi=hindi+1;

    [Phivvvvc_vec,model.hind{hindi}]=chain4c_tensor(Phiz_vec,Phizz_vec,Phizzz_vec,Phizzzz_vec,...
        PhiPIv_vec,PhiPIvvc_vec,PhiPIvvvc_vec,PhiPIvvvvc_vec,...
        model.Phivars_chain4c_M2,model.Phivars_chain4c_M3,model.Phivars_chain4c_M4,...
        model.hind{hindi},n_ind,maxload,'vec',model.tilPhiz,model.maxtilPhiz);
    hindi=hindi+1;
    
    [Phivvvvvc_vec,model.hind4{hindi4}]=chain5c_tensor(Phiz_vec,Phizz_vec,Phizzz_vec,Phizzzz_vec,Phizzzzz_vec,...
        PhiPIv_vec,PhiPIvvc_vec,PhiPIvvvc_vec,PhiPIvvvvc_vec,PhiPIvvvvvc_vec,...
        model.Phivars_chain5c_M1,model.Phivars_chain5c_M2,model.Phivars_chain5c_M3,model.Phivars_chain5c_M4,model.Phivars_chain5c_M5,model.Phivars_chain5c_M6,...
        model.hind4{hindi4},n_ind,maxload,'vec',model.tilPhiz,model.maxtilPhiz);
    hindi4=hindi4+1;

%     stochfv_vec.vals=multcol(stochfv_vec.vals,Pvec);
%     stochfvvc_vec.vals=multcol(stochfvvc_vec.vals,Pvec);
%     stochfvvvc_vec.vals=multcol(stochfvvvc_vec.vals,Pvec);
%     stochfvvvvc_vec.vals=multcol(stochfvvvvc_vec.vals,Pvec);
%     stochfvvvvvc_vec.vals=multcol(stochfvvvvvc_vec.vals,Pvec);

    [ Phivv_vec,model.hind{hindi} ] = uncompressderivs( Phivvc_vec,2,n_Phivars,model.Phiv(:,Phivars),model.hind{hindi} );
    hindi=hindi+1;
    clear Phivvc_vec

    [ Phivvv_vec,model.hind{hindi} ] = uncompressderivs( Phivvvc_vec,3,n_Phivars,model.Phiv(:,Phivars),model.hind{hindi} );
    hindi=hindi+1;
    clear Phivvvc_vec

    [ Phivvvv_vec,model.hind{hindi} ] = uncompressderivs( Phivvvvc_vec,4,n_Phivars,model.Phiv(:,Phivars),model.hind{hindi} );
    hindi=hindi+1;
    clear Phivvvvc_vec

    [ Phivvvvv_vec,model.hind4{hindi4} ] = uncompressderivs( Phivvvvvc_vec,5,n_Phivars,model.Phiv(:,Phivars),model.hind4{hindi4} );
    hindi4=hindi4+1;
    clear Phivvvvvc_vec

end
tothindi=tothindi+8+count_Phi_n*(1+3+4+6)+1+4+6+9;
hindi=tothindi;