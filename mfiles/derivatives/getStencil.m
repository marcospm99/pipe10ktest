function [fStencil,nfStencil] = getStencil(n,nder)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Funci�n que devuelve la matriz de nodos de influencia relativa al nodo
% central del stencil
% Alfonso Pallares - 23-Sept-2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stencil f

if mod(n,2)==0
    error('n must be odd')
end

if mod(nder,2)==0
    error('n must be odd')
end


f = ones(1,n);
%--
[~,n] = size(f);
n1    = (n+1)/2;
[~,b] = find(f==1);
fStencil = b-n1;

%% Stencil f^{(n)}
nf = ones(1,nder);
%--
[~,m] = size(nf);
m1    = (m+1)/2;
[~,a] = find(nf==1);
nfStencil = a-m1;