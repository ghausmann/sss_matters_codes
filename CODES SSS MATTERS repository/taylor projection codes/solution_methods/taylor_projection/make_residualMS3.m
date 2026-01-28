function make_residualMS3(f,Phi,yp,y,xp,x,shocks,chip,chi,symparams,order,ufun,u,uu0,logicalparams,logicalvars)
% make_residualMS(f,Phi,yp,y,xp,x,symparams_msp,symparams_ms,symparams,eta,order,auxfuns,auxvars,model.uu0,logicalparams,logicalvars);
% The function creates a Matlab file that computes the residual function for
% a model with Markov-switching parameters.

all_symparams=[symparams(:);chi(:);chip(:)];

make_residual3(f,Phi,yp,y,xp,x,shocks,all_symparams,order,ufun,u,uu0,'MS',logicalparams,logicalvars);




