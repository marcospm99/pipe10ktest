function [nfDistance] = nfStencilPivot(nfDistance0,i)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Esta funci�n elimina la fila de informaci�n que corresponde al nodo de
% estudio. Se elimina en base al pivotaje potserior que sobre \alpha_{i} se
% realiza en la resoluci�n del SEL
% Alfonso Pallar�s - 23-Sept-2014
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