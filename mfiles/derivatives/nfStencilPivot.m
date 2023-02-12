function [nfDistance] = nfStencilPivot(nfDistance0,i)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Esta función elimina la fila de información que corresponde al nodo de
% estudio. Se elimina en base al pivotaje potserior que sobre \alpha_{i} se
% realiza en la resolución del SEL
% Alfonso Pallarés - 23-Sept-2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
A         = nfDistance0(:,1)-nfDistance0(:,2);
rowDelete = find(A==0);
if isempty(rowDelete)
    error('Function: nfStencilPivot(nfDistance0,i) -> Error: Study node not found in matrix nfDistance. Whole code check is required. If this is correct: re-build current and previous functions');
end
nfDistance0(rowDelete,:) = [];
%--
nfDistance = nfDistance0;