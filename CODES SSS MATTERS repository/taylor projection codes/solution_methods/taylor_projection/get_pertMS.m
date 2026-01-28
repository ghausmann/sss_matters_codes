function [derivs,stoch_pert,nonstoch_pert,model]=get_pertMS(model,params,msparams_ss,M,nxss,nyss,algo)
%Get a standard perturbation solution for fixed parameters.

[derivs,stoch_pert,nonstoch_pert,model]=get_pert(model,[params(:);msparams_ss(:);msparams_ss(:)],M,nxss,nyss,algo);

end

