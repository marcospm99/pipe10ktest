function u = dirichletCC(A,B,f,CC1,CC2,nx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Función para calcular la solución numérica aplicando condiciones de
% contorno de Dirichlet
% Alfonso Pallarés - 23-Sept-2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construcción del vector $f$
u = [CC1;...
    zeros(nx-2,1);...
    CC2];
% Resolución del SEl
g          = A*f - B*u;
u(2:end-1) = B(2:end-1,2:end-1)\g(2:end-1);
