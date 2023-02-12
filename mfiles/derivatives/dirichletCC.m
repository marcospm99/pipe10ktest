function u = dirichletCC(A,B,f,CC1,CC2,nx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Funci�n para calcular la soluci�n num�rica aplicando condiciones de
% contorno de Dirichlet
% Alfonso Pallar�s - 23-Sept-2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construcci�n del vector $f$
u = [CC1;...
    zeros(nx-2,1);...
    CC2];
% Resoluci�n del SEl
g          = A*f - B*u;
u(2:end-1) = B(2:end-1,2:end-1)\g(2:end-1);
