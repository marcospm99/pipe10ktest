function [fStencil,nfStencil] = getStencilX()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Función que devuelve la matriz de nodos de influencia relativa al nodo
% central del stencil
% Alfonso Pallares - 23-Sept-2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stencil f
f = [1,1,1,1,1];
%--
[~,n] = size(f);
n1    = (n+1)/2;
[~,b] = find(f==1);
fStencil = b-n1;

%% Stencil f^{(n)}
nf = [1,1,1];
%--
[~,m] = size(nf);
m1    = (m+1)/2;
[~,a] = find(nf==1);
nfStencil = a-m1;