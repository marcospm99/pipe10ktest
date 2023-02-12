function [fMat,nfMat,lf,lnf] = taylorExpansion(fDistance,nfDistance,derivada)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Funci�n que devuelve el desarrollo de Taylor unidimensional ajustando el
% orden del desarrollo a la informaci�n disponible, es decir, funci�n del
% stencil de c�lculo.
% Alfonso Pallar�s - 23-Sept-2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preparaci�n del desarrollo
hX    = transpose(fDistance(:,3));
hXX   = transpose(nfDistance(:,3));
const = derivada+1;
%
lf    = length(hX);
lnf   = length(hXX);
fMat  = [ones(1,lf);...
         zeros(lf+lnf-1,lf)];
nfMat = [zeros(derivada,lnf);...
         -ones(1,lnf);...
         zeros(lf+lnf-const,lnf)];
%% Desarrollo de los t�rminos de la funci�n
for k = 1:(lf+lnf-1)
    fMat(k+1,:) = hX.^(k)/factorial(k);
end
%% Desarrollo de los t�rminos de la SEGUNDA derivada de la funci�n
for k = 1:(lf+lnf-const)
    nfMat(k+const,:) = -hXX.^(k)/factorial(k);
end