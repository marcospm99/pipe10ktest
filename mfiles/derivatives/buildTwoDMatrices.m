function [A,B,C,D,fx,fppx,fy,fppy] = buildTwoDMatrices(x,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Función que devuelve las matrices del esquema CFD asociadas a cada eje
% del acoplamiento.
% Alfonso Pallarés - 23-Sept-2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Los vectores x e y entran como dato a la función
derivada = 2;
%% Construcción de las matrices A y B
[fxStencil,nfxStencil] = getStencilX();
[A,B] = matrixCFD(x,fxStencil,nfxStencil,derivada);
fx    = 1+2*fxStencil(1,end);
fppx  = 1+2*nfxStencil(1,end);
%% Construcción de las matrices A y B
[fyStencil,nfyStencil] = getStencilY();
[C,D] = matrixCFD(y,fyStencil,nfyStencil,derivada);
fy    = 1+2*fyStencil(1,end);
fppy  = 1+2*nfyStencil(1,end);