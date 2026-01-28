function [n,uu0]=find_n_new(pi_,u)
% Given u=pi_(v,u), how many compositions of type pi_(v,pi_(v,u)) are needed to get a
% function of v only.
%
% © Copyright, Oren Levintal, September 12, 2017.

if isempty(pi_)
    n=0;
    uu0=zeros(0,0);
else
    n_u=length(u);
    if length(pi_)~=n_u
        error('wrong size of function')
    end
    pi_=pi_(:);
    uu0=double(logical(jacobian(pi_,u)~=0));
    if max(abs(real(eig(uu0))))>=1-2*eps
        error('Auxiliary functions are not defined properly.')
    end
    n=1;
    uu=uu0;
    nuu=uu*uu0;
    while ~isequal(logical(nuu),logical(uu))
        n=n+1;
        uu=nuu;
        nuu=uu*uu0;
    end
end

end