function make_residualMS(f,Phi,yp,y,xp,x,chip,chi,symparams,eta,order,ufun,u,uu0,logicalparams,logicalvars)
% make_residualMS(f,Phi,yp,y,xp,x,symparams_msp,symparams_ms,symparams,eta,order,auxfuns,auxvars,model.uu0,logicalparams,logicalvars);
% The function creates a Matlab file that computes the residual function for
% a model with Markov-switching parameters.

all_symparams=[symparams(:);chi(:);chip(:)];

make_residual(f,Phi,yp,y,xp,x,all_symparams,eta,order,ufun,u,uu0,'MS',logicalparams,logicalvars);




