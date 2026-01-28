if ~isempty(model.prefrows)

    if approx>=1
        prefv=changecols(prefv,model.prefvars,n_v_tp,1);
        tempv=1:n_v;
        tempv(2*n_y+n_x)=[];
        tempv(end)=[];
        prefv=changecols(prefv,tempv,n_v,1);
        prefv=changerows(prefv,model.prefrows,n_f);

        stochfv_vec=changecols(stochfv_vec,model.stochfvars,n_v_tp,1);
        stochfv_vec=changecols(stochfv_vec,tempv,n_v,1);
        stochfv_vec=changerows(stochfv_vec,model.stochfrows,n_f);

        Phiv_vec=changecols(Phiv_vec,model.Phivars,n_v_tp,1);
        Phiv_vec=changecols(Phiv_vec,tempv,n_v,1);
        
        [prei,prej,prevals]=tfind(unfold(prefv));
        [stochi,stochj,stochvals]=tfind(unfold(stochfv_vec));
        [Phii,Phij,Phivals]=tfind(unfold(Phiv_vec));
        Phii=Phii+n_f-n_x2;

        fv=sparse(double([prei;stochi;Phii]),double([prej;stochj;Phij]),double([prevals(:);stochvals(:);Phivals(:)]),n_f,n_v);
    end
    if approx>=2
        prefvv=changecols(prefvv,model.prefvars,n_v_tp,1);
        prefvv=changecols(prefvv,model.prefvars,n_v_tp,2);
        prefvv=changecols(prefvv,tempv,n_v,1);
        prefvv=changecols(prefvv,tempv,n_v,2);
        prefvv=changerows(prefvv,model.prefrows,n_f);

        stochfvv_vec=changecols(stochfvv_vec,model.stochfvars,n_v_tp,1);
        stochfvv_vec=changecols(stochfvv_vec,model.stochfvars,n_v_tp,2);
        stochfvv_vec=changecols(stochfvv_vec,tempv,n_v,1);
        stochfvv_vec=changecols(stochfvv_vec,tempv,n_v,2);
        stochfvv_vec=changerows(stochfvv_vec,model.stochfrows,n_f);

        Phivv_vec=changecols(Phivv_vec,model.Phivars,n_v_tp,1);
        Phivv_vec=changecols(Phivv_vec,model.Phivars,n_v_tp,2);
        Phivv_vec=changecols(Phivv_vec,tempv,n_v,1);
        Phivv_vec=changecols(Phivv_vec,tempv,n_v,2);
        
        [prei,prej,prevals]=tfind(unfold(prefvv));
        [stochi,stochj,stochvals]=tfind(unfold(stochfvv_vec));
        [Phii,Phij,Phivals]=tfind(unfold(Phivv_vec));
        Phii=Phii+n_f-n_x2;

        alli=double([prei;stochi;Phii]);
        allj=double([prej;stochj;Phij]);
        allvals=double([prevals(:);stochvals(:);Phivals(:)]);

        newi=sub2ind([n_f,n_v^2],alli,allj);
        newj=ones(size(newi));

        fvv=sparse(newi,newj,allvals,n_f*n_v^2,1);
    end
    if approx>=3
        prefvvv=changecols(prefvvv,model.prefvars,n_v_tp,1);
        prefvvv=changecols(prefvvv,model.prefvars,n_v_tp,2);
        prefvvv=changecols(prefvvv,model.prefvars,n_v_tp,3);
        prefvvv=changecols(prefvvv,tempv,n_v,1);
        prefvvv=changecols(prefvvv,tempv,n_v,2);
        prefvvv=changecols(prefvvv,tempv,n_v,3);
        prefvvv=changerows(prefvvv,model.prefrows,n_f);

        stochfvvv_vec=changecols(stochfvvv_vec,model.stochfvars,n_v_tp,1);
        stochfvvv_vec=changecols(stochfvvv_vec,model.stochfvars,n_v_tp,2);
        stochfvvv_vec=changecols(stochfvvv_vec,model.stochfvars,n_v_tp,3);
        stochfvvv_vec=changecols(stochfvvv_vec,tempv,n_v,1);
        stochfvvv_vec=changecols(stochfvvv_vec,tempv,n_v,2);
        stochfvvv_vec=changecols(stochfvvv_vec,tempv,n_v,3);
        stochfvvv_vec=changerows(stochfvvv_vec,model.stochfrows,n_f);

        Phivvv_vec=changecols(Phivvv_vec,model.Phivars,n_v_tp,1);
        Phivvv_vec=changecols(Phivvv_vec,model.Phivars,n_v_tp,2);
        Phivvv_vec=changecols(Phivvv_vec,model.Phivars,n_v_tp,3);
        Phivvv_vec=changecols(Phivvv_vec,tempv,n_v,1);
        Phivvv_vec=changecols(Phivvv_vec,tempv,n_v,2);
        Phivvv_vec=changecols(Phivvv_vec,tempv,n_v,3);

        [prei,prej,prevals]=tfind(unfold(prefvvv));
        [stochi,stochj,stochvals]=tfind(unfold(stochfvvv_vec));
        [Phii,Phij,Phivals]=tfind(unfold(Phivvv_vec));
        Phii=Phii+n_f-n_x2;

        alli=double([prei;stochi;Phii]);
        allj=double([prej;stochj;Phij]);
        allvals=double([prevals(:);stochvals(:);Phivals(:)]);

        newi=sub2ind([n_f,n_v^3],alli,allj);
        newj=ones(size(newi));

        fvvv=sparse(newi,newj,allvals,n_f*n_v^3,1);
    end
    if approx>=4
        prefvvvv=changecols(prefvvvv,model.prefvars,n_v_tp,1);
        prefvvvv=changecols(prefvvvv,model.prefvars,n_v_tp,2);
        prefvvvv=changecols(prefvvvv,model.prefvars,n_v_tp,3);
        prefvvvv=changecols(prefvvvv,model.prefvars,n_v_tp,4);
        prefvvvv=changecols(prefvvvv,tempv,n_v,1);
        prefvvvv=changecols(prefvvvv,tempv,n_v,2);
        prefvvvv=changecols(prefvvvv,tempv,n_v,3);
        prefvvvv=changecols(prefvvvv,tempv,n_v,4);
        prefvvvv=changerows(prefvvvv,model.prefrows,n_f);

        stochfvvvv_vec=changecols(stochfvvvv_vec,model.stochfvars,n_v_tp,1);
        stochfvvvv_vec=changecols(stochfvvvv_vec,model.stochfvars,n_v_tp,2);
        stochfvvvv_vec=changecols(stochfvvvv_vec,model.stochfvars,n_v_tp,3);
        stochfvvvv_vec=changecols(stochfvvvv_vec,model.stochfvars,n_v_tp,4);
        stochfvvvv_vec=changecols(stochfvvvv_vec,tempv,n_v,1);
        stochfvvvv_vec=changecols(stochfvvvv_vec,tempv,n_v,2);
        stochfvvvv_vec=changecols(stochfvvvv_vec,tempv,n_v,3);
        stochfvvvv_vec=changecols(stochfvvvv_vec,tempv,n_v,4);
        stochfvvvv_vec=changerows(stochfvvvv_vec,model.stochfrows,n_f);

        Phivvvv_vec=changecols(Phivvvv_vec,model.Phivars,n_v_tp,1);
        Phivvvv_vec=changecols(Phivvvv_vec,model.Phivars,n_v_tp,2);
        Phivvvv_vec=changecols(Phivvvv_vec,model.Phivars,n_v_tp,3);
        Phivvvv_vec=changecols(Phivvvv_vec,model.Phivars,n_v_tp,4);
        Phivvvv_vec=changecols(Phivvvv_vec,tempv,n_v,1);
        Phivvvv_vec=changecols(Phivvvv_vec,tempv,n_v,2);
        Phivvvv_vec=changecols(Phivvvv_vec,tempv,n_v,3);
        Phivvvv_vec=changecols(Phivvvv_vec,tempv,n_v,4);

        [prei,prej,prevals]=tfind(unfold(prefvvvv));
        [stochi,stochj,stochvals]=tfind(unfold(stochfvvvv_vec));
        [Phii,Phij,Phivals]=tfind(unfold(Phivvvv_vec));
        Phii=Phii+n_f-n_x2;

        alli=double([prei;stochi;Phii]);
        allj=double([prej;stochj;Phij]);
        allvals=double([prevals(:);stochvals(:);Phivals(:)]);

        newi=sub2ind([n_f,n_v^4],alli,allj);
        newj=ones(size(newi));

        fvvvv=sparse(newi,newj,allvals,n_f*n_v^4,1);
    end
else
    if approx>=1
%         prefv=changecols(prefv,model.prefvars,n_v_tp,1);
        tempv=1:n_v;
        tempv(2*n_y+n_x)=[];
        tempv(end)=[];
%         prefv=changecols(prefv,tempv,n_v,1);
%         prefv=changerows(prefv,model.prefrows,n_f);

        stochfv_vec=changecols(stochfv_vec,model.stochfvars,n_v_tp,1);
        stochfv_vec=changecols(stochfv_vec,tempv,n_v,1);
        stochfv_vec=changerows(stochfv_vec,model.stochfrows,n_f);

        Phiv_vec=changecols(Phiv_vec,model.Phivars,n_v_tp,1);
        Phiv_vec=changecols(Phiv_vec,tempv,n_v,1);

%         [prei,prej,prevals]=tfind(unfold(prefv));
        [stochi,stochj,stochvals]=tfind(unfold(stochfv_vec));
        [Phii,Phij,Phivals]=tfind(unfold(Phiv_vec));
        Phii=Phii+n_f-n_x2;

%         fv=sparse(double([prei;stochi]),double([prej;stochj]),double([prevals(:);stochvals(:)]),n_f,n_v);
        fv=sparse(double([stochi;Phii]),double([stochj;Phij]),double([stochvals(:);Phivals(:)]),n_f,n_v);
    end
    if approx>=2
%         prefvv=changecols(prefvv,model.prefvars,n_v_tp,1);
%         prefvv=changecols(prefvv,model.prefvars,n_v_tp,2);
%         prefvv=changecols(prefvv,tempv,n_v,1);
%         prefvv=changecols(prefvv,tempv,n_v,2);
%         prefvv=changerows(prefvv,model.prefrows,n_f);

        stochfvv_vec=changecols(stochfvv_vec,model.stochfvars,n_v_tp,1);
        stochfvv_vec=changecols(stochfvv_vec,model.stochfvars,n_v_tp,2);
        stochfvv_vec=changecols(stochfvv_vec,tempv,n_v,1);
        stochfvv_vec=changecols(stochfvv_vec,tempv,n_v,2);
        stochfvv_vec=changerows(stochfvv_vec,model.stochfrows,n_f);

        Phivv_vec=changecols(Phivv_vec,model.Phivars,n_v_tp,1);
        Phivv_vec=changecols(Phivv_vec,model.Phivars,n_v_tp,2);
        Phivv_vec=changecols(Phivv_vec,tempv,n_v,1);
        Phivv_vec=changecols(Phivv_vec,tempv,n_v,2);
        

%         [prei,prej,prevals]=tfind(unfold(prefvv));
        [stochi,stochj,stochvals]=tfind(unfold(stochfvv_vec));
        [Phii,Phij,Phivals]=tfind(unfold(Phivv_vec));
        Phii=Phii+n_f-n_x2;

        alli=double([stochi;Phii]);
        allj=double([stochj;Phij]);
        allvals=double([stochvals(:);Phivals(:)]);

        newi=sub2ind([n_f,n_v^2],alli,allj);
        newj=ones(size(newi));

        fvv=sparse(newi,newj,allvals,n_f*n_v^2,1);
    end
    if approx>=3
%         prefvvv=changecols(prefvvv,model.prefvars,n_v_tp,1);
%         prefvvv=changecols(prefvvv,model.prefvars,n_v_tp,2);
%         prefvvv=changecols(prefvvv,model.prefvars,n_v_tp,3);
%         prefvvv=changecols(prefvvv,tempv,n_v,1);
%         prefvvv=changecols(prefvvv,tempv,n_v,2);
%         prefvvv=changecols(prefvvv,tempv,n_v,3);
%         prefvvv=changerows(prefvvv,model.prefrows,n_f);

        stochfvvv_vec=changecols(stochfvvv_vec,model.stochfvars,n_v_tp,1);
        stochfvvv_vec=changecols(stochfvvv_vec,model.stochfvars,n_v_tp,2);
        stochfvvv_vec=changecols(stochfvvv_vec,model.stochfvars,n_v_tp,3);
        stochfvvv_vec=changecols(stochfvvv_vec,tempv,n_v,1);
        stochfvvv_vec=changecols(stochfvvv_vec,tempv,n_v,2);
        stochfvvv_vec=changecols(stochfvvv_vec,tempv,n_v,3);
        stochfvvv_vec=changerows(stochfvvv_vec,model.stochfrows,n_f);

        Phivvv_vec=changecols(Phivvv_vec,model.Phivars,n_v_tp,1);
        Phivvv_vec=changecols(Phivvv_vec,model.Phivars,n_v_tp,2);
        Phivvv_vec=changecols(Phivvv_vec,model.Phivars,n_v_tp,3);
        Phivvv_vec=changecols(Phivvv_vec,tempv,n_v,1);
        Phivvv_vec=changecols(Phivvv_vec,tempv,n_v,2);
        Phivvv_vec=changecols(Phivvv_vec,tempv,n_v,3);


%         [prei,prej,prevals]=tfind(unfold(prefvvv));
        [stochi,stochj,stochvals]=tfind(unfold(stochfvvv_vec));
        [Phii,Phij,Phivals]=tfind(unfold(Phivvv_vec));
        Phii=Phii+n_f-n_x2;

        alli=double([stochi;Phii]);
        allj=double([stochj;Phij]);
        allvals=double([stochvals(:);Phivals(:)]);

        newi=sub2ind([n_f,n_v^3],alli,allj);
        newj=ones(size(newi));

        fvvv=sparse(newi,newj,allvals,n_f*n_v^3,1);
    end
    if approx>=4
%         prefvvvv=changecols(prefvvvv,model.prefvars,n_v_tp,1);
%         prefvvvv=changecols(prefvvvv,model.prefvars,n_v_tp,2);
%         prefvvvv=changecols(prefvvvv,model.prefvars,n_v_tp,3);
%         prefvvvv=changecols(prefvvvv,model.prefvars,n_v_tp,4);
%         prefvvvv=changecols(prefvvvv,tempv,n_v,1);
%         prefvvvv=changecols(prefvvvv,tempv,n_v,2);
%         prefvvvv=changecols(prefvvvv,tempv,n_v,3);
%         prefvvvv=changecols(prefvvvv,tempv,n_v,4);
%         prefvvvv=changerows(prefvvvv,model.prefrows,n_f);

        stochfvvvv_vec=changecols(stochfvvvv_vec,model.stochfvars,n_v_tp,1);
        stochfvvvv_vec=changecols(stochfvvvv_vec,model.stochfvars,n_v_tp,2);
        stochfvvvv_vec=changecols(stochfvvvv_vec,model.stochfvars,n_v_tp,3);
        stochfvvvv_vec=changecols(stochfvvvv_vec,model.stochfvars,n_v_tp,4);
        stochfvvvv_vec=changecols(stochfvvvv_vec,tempv,n_v,1);
        stochfvvvv_vec=changecols(stochfvvvv_vec,tempv,n_v,2);
        stochfvvvv_vec=changecols(stochfvvvv_vec,tempv,n_v,3);
        stochfvvvv_vec=changecols(stochfvvvv_vec,tempv,n_v,4);
        stochfvvvv_vec=changerows(stochfvvvv_vec,model.stochfrows,n_f);

        Phivvvv_vec=changecols(Phivvvv_vec,model.Phivars,n_v_tp,1);
        Phivvvv_vec=changecols(Phivvvv_vec,model.Phivars,n_v_tp,2);
        Phivvvv_vec=changecols(Phivvvv_vec,model.Phivars,n_v_tp,3);
        Phivvvv_vec=changecols(Phivvvv_vec,model.Phivars,n_v_tp,4);
        Phivvvv_vec=changecols(Phivvvv_vec,tempv,n_v,1);
        Phivvvv_vec=changecols(Phivvvv_vec,tempv,n_v,2);
        Phivvvv_vec=changecols(Phivvvv_vec,tempv,n_v,3);
        Phivvvv_vec=changecols(Phivvvv_vec,tempv,n_v,4);


%         [prei,prej,prevals]=tfind(unfold(prefvvvv));
        [stochi,stochj,stochvals]=tfind(unfold(stochfvvvv_vec));
        [Phii,Phij,Phivals]=tfind(unfold(Phivvvv_vec));
        Phii=Phii+n_f-n_x2;

        alli=double([stochi;Phii]);
        allj=double([stochj;Phij]);
        allvals=double([stochvals(:);Phivals(:)]);

        newi=sub2ind([n_f,n_v^4],alli,allj);
        newj=ones(size(newi));

        fvvvv=sparse(newi,newj,allvals,n_f*n_v^4,1);
    end    
end
