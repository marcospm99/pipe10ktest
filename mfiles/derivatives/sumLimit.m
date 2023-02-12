function[limit] = sumLimit(h,nodes,M)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%       27-Nov-2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Esta funci�n devuelve los l�mites que se tienen que imponer en el
% sumatorio 
%
% h -> impuesto desde fuera. Es la fila de la matriz resultado
% nodes -> n�mero de nodos que se emplean para calcular M
% M -> tama�o de la matriz que impone la restricci�n
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nS  = (nodes-1)/2;
if h<(nS+1)
    limit = [1,h+nS];
elseif (h>nS)&&(h<(M-nS))
    limit = [h-nS,h+nS];
elseif h>(M-nS-1)
    limit = [h-nS,M];
end

