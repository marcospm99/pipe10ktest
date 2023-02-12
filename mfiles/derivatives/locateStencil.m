function [fStencilMesh,nfStencilMesh] = locateStencil(nodo,nx,fStencil,nfStencil)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Esta funci�n permite colocar correctamente el stencil definido dentro de
% la malla. Localiza los puntos del stencil que, al acercarse a la pared,
% se quedan fuera y los elimina para que �nicmaente prevalezcan los nodos
% del stencil que quedan estrictamente dentro de la malla.
%
% Esta funci�n precisa como entrada:
%  - nodo      -> N�mero del nodo de estudio
%  - nx        -> Tama�o de la malla
%  - fStencil  -> Stencil de la funci�n   
%  - nfStencil -> Stencil de la n-�sima derivada
% Esta funci�n devuelve como salida:
%  - fStencilMesh  -> Matriz con la informaci�n del stencil en la malla
%  - nfStencilMesh -> Matriz con la informaci�n del stencil en la malla
%
% Alfonso Pallar�s - 23-Sept-2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Primera localizaci�n del stencil dentro de la malla
fStencilMesh  = fStencil+nodo;
nfStencilMesh = nfStencil+nodo;
% Borrado de los nodos con componente <= (0) � > (nx)
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
