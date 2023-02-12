function [fDistance,nfDistance,s1,s2] = nodalDistance(fStencilMesh,nfStencilMesh,x,px)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Función que devuelve las distancias entre los nodos relacionados a través
% del stencil
% Alfonso Pallarés - 23-Sept-2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fStencil
[~,s1]    = size(fStencilMesh);
fDistance = zeros(s1,3);%[x0,x1,f_hx]->h_i=i_1-i_0
for i=1:s1
    x1f = fStencilMesh(1,i);
    f_hx  = x(x1f)-x(px);
    fDistance(i,:) = [px,x1f,f_hx];
end
% nfStencil
[~,s2]     = size(nfStencilMesh);
nfDistance = zeros(s2,3);%[x0,x1,f_hx]->h_i=i_1-i_0
for i=1:s2
    x1nf   = nfStencilMesh(1,i);
    nf_hx  = x(x1nf)-x(px);
    nfDistance(i,:) = [px,x1nf,nf_hx];
end