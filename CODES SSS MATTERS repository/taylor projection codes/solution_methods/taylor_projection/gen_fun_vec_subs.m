function gen_fun_vec_subs(f,symparams,v,ufun,u,fname,varargin)
%
% © Copyright, Oren Levintal, September 11, 2017.

rowformat=0; % function returns a column vector (number of columns>1 in vectorized form)
if ~isempty(varargin)
    if strcmp(varargin{1},'row')
        rowformat=1; % function returns a row vector (number of rows>1 in vectorized form)
    end
end
if nargin==8
    uu0=varargin{2};
end

f=f(:);
if numel(f)>1
    nnzf=1-logical(f==0);
    index=find(nnzf);
else
    index=1;
end
fid = fopen([fname '_fun.m'], 'w');
fprintf(fid,'%s\n', ['function function_f=' fname '_fun(variables_v,parameters)']);
for i=1:length(symparams)
    if strcmp(char(symparams(i)),'function_f') || strcmp(char(symparams(i)),'variables_v') || strcmp(char(symparams(i)),'parameters')
        error([char(symparams(i)) 'is reserved. Change variable/parameter name.']);
    end
    if rowformat==0
        fprintf(fid,'%s\n', [char(symparams(i)) '=parameters(' num2str(i) ',:);']);
    elseif rowformat==1
        fprintf(fid,'%s\n', [char(symparams(i)) '=parameters(:,' num2str(i) ');']);
    end
end
for i=1:length(v)
    if strcmp(char(v(i)),'function_f') || strcmp(char(v(i)),'variables_v') || strcmp(char(v(i)),'parameters')
        error([char(v(i)) 'is reserved. Change variable/parameter name.']);
    end
    if rowformat==0
        fprintf(fid,'%s\n', [char(v(i)) '=variables_v(' num2str(i) ',:);']);
    elseif rowformat==1
        fprintf(fid,'%s\n', [char(v(i)) '=variables_v(:,' num2str(i) ');']);
    end
end

%%%%%%% substituted vars - version 1

% subsu=ufun;
% newsubsu=subs(subsu,u,ufun);
% 
% n_u=length(u);
% A=cell(n_u,1);
% l=1;
% uind=1:n_u;
% while ~isequal(subsu,newsubsu)
%     for i=uind(logical(newsubsu==subsu))
%         texp=char(ufun(i));
%         texp=strrep(texp,'*','.*');
%         texp=strrep(texp,'/','./');
%         texp=strrep(texp,'^','.^');
%         A{l}=[char(u(i)) '=' texp ';'];l=l+1;
%     end
%     uind=uind(logical(newsubsu~=subsu));
%     subsu=newsubsu(logical(newsubsu~=subsu));
%     newsubsu=subs(subsu,u,ufun);
% end
% for i=uind
%     texp=char(ufun(i));
%     texp=strrep(texp,'*','.*');
%     texp=strrep(texp,'/','./');
%     texp=strrep(texp,'^','.^');
%     A{l}=[char(u(i)) '=' texp ';'];l=l+1;
% end
% 
% for i=1:l-1
%     fprintf(fid,'%s\n', A{i});
% end

%%%%%%% substituted vars - version 2

if nargin~=8
    uu0=double(logical(jacobian(ufun(:),u)~=0));
end

uu=uu0;

n_u=length(u);
A=cell(n_u,1);
l=1;
uind=1:n_u;
while nnz(uu)~=0
    for i=uind(logical(sum(uu,2)==0))
        texp=char(ufun(i));
        texp=strrep(texp,'*','.*');
        texp=strrep(texp,'/','./');
        texp=strrep(texp,'^','.^');
        A{l}=[char(u(i)) '=' texp ';'];l=l+1;
    end
    uind=uind(logical(sum(uu,2)~=0));
    uu=uu(logical(sum(uu,2)~=0),:);
    uu=uu*uu0;
end
for i=uind
    texp=char(ufun(i));
    texp=strrep(texp,'*','.*');
    texp=strrep(texp,'/','./');
    texp=strrep(texp,'^','.^');
    A{l}=[char(u(i)) '=' texp ';'];l=l+1;
end

for i=1:l-1
    fprintf(fid,'%s\n', A{i});
end

%%%%%%%%%%%


if isempty(f)
    if rowformat==0
        fprintf(fid,'%s\n', ['function_f=zeros(0,size(variables_v,2));']);
    elseif rowformat==1
        fprintf(fid,'%s\n', ['function_f=zeros(size(variables_v,1),0);']);
    end
else
    if rowformat==0
        fprintf(fid,'%s\n', ['function_f=zeros(' num2str(numel(f)) ',size(variables_v,2));']);
    elseif rowformat==1
        fprintf(fid,'%s\n', ['function_f=zeros(size(variables_v,1),' num2str(numel(f)) ');']);
    end
    for i=index(:)'
        if rowformat==0
            disp_fun_vec('function_f',f(i),i,fid);
        elseif rowformat==1
            disp_fun_vec_row('function_f',f(i),i,fid);
        end
    end
end
fclose(fid);
