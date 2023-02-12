function [A,b] = buildMatrix(fMat,nfMat,derivada)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Función para la construcción y condicionamiento de la matriz para el
% cálculo de los coeficientes del esquema de Diferencias Finitas Compactas
% Alfonso Pallarés - 23-Sept-2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construcción de la matriz y vector
AnoCond = [nfMat,fMat];
[m,n]   = size(AnoCond);
const   = derivada+1;
if m~=n
    error('Function: buildMatrix(fMat,nfMat) -> Error: Resultant matrix ([fMat,nfMat]) NOT square. Check previous functions');
end
bnoCond = [zeros(derivada,1);...
           ones(1,1);...
           zeros(m-const,1)];
A = zeros(size(AnoCond));
b = zeros(size(bnoCond));
%% Condicionamiento
for i = 1:m
    vec1 = AnoCond(i,:);
    vec2 = sort(abs(vec1));
    vec3 = vec2(vec2~=0);
    cond = vec3(1);
    A(i,:) = AnoCond(i,:)/cond;
    b(i,1) = bnoCond(i,1)/cond;
end