function [x] = getDiscretization(a,b,n,flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Funci�n que devuelve el vector del campo discreto \vec{x}. La
% discretizaci�n se realiza sobre el campo continuo [a,b] \in\mathbb{R}^{1}
% y contiene n puntos.
% Es necesario indicar la ley que se emplea para discretizar
% Alfonso Pallares - 23-Sept-2014
% flag = 1 --> equiespaciado
% Flag = 2 --> no equiespaciado
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if b<a
    error('Function: getDiscretization(a,b,n) -> Error: b must be greater than a')
end

if flag==1
    x = linspace(a,b,n)';
else
    x = vectorNoEquiespaciado(n);
end