function [fMat,nfMat,lf,lnf] = taylorExpansion(fDistance,nfDistance,derivada)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Función que devuelve el desarrollo de Taylor unidimensional ajustando el
% orden del desarrollo a la información disponible, es decir, función del
% stencil de cálculo.
% Alfonso Pallarés - 23-Sept-2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preparación del desarrollo
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
%% Desarrollo de los términos de la función
for k = 1:(lf+lnf-1)
    fMat(k+1,:) = hX.^(k)/factorial(k);
end
%% Desarrollo de los términos de la SEGUNDA derivada de la función
for k = 1:(lf+lnf-const)
    nfMat(k+const,:) = -hXX.^(k)/factorial(k);
end