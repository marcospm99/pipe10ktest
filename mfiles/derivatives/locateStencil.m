function [fStencilMesh,nfStencilMesh] = locateStencil(nodo,nx,fStencil,nfStencil)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Esta función permite colocar correctamente el stencil definido dentro de
% la malla. Localiza los puntos del stencil que, al acercarse a la pared,
% se quedan fuera y los elimina para que únicmaente prevalezcan los nodos
% del stencil que quedan estrictamente dentro de la malla.
%
% Esta función precisa como entrada:
%  - nodo      -> Número del nodo de estudio
%  - nx        -> Tamaño de la malla
%  - fStencil  -> Stencil de la función   
%  - nfStencil -> Stencil de la n-ésima derivada
% Esta función devuelve como salida:
%  - fStencilMesh  -> Matriz con la información del stencil en la malla
%  - nfStencilMesh -> Matriz con la información del stencil en la malla
%
% Alfonso Pallarés - 23-Sept-2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Primera localización del stencil dentro de la malla
fStencilMesh  = fStencil+nodo;
nfStencilMesh = nfStencil+nodo;
% Borrado de los nodos con componente <= (0) ó > (nx)
% fStencil
[~,s2f] = find(fStencilMesh<=0);
s2fp      = unique(transpose(s2f));
fStencilMesh(:,s2fp) = [];
[~,s2f] = find(fStencilMesh(1,:)>nx);
fStencilMesh(:,transpose(s2f)) = [];
% nfStencil
[~,s2nf] = find(nfStencilMesh<=0);
s2nfp       = unique(transpose(s2nf));
nfStencilMesh(:,s2nfp) = [];
[~,s2nf] = find(nfStencilMesh(1,:)>nx);
nfStencilMesh(:,transpose(s2nf)) = [];
