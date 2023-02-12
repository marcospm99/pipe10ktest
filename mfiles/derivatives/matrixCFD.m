function [A,B,Ad,Bd] = matrixCFD(x,fStencil,nfStencil,derivada,flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Funci�n para el c�lculo de las matrices A y B del esquema de diferencias
% finitas compactas
% Alfonso Pallares - 23-Sept-2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nx = length(x);
B = zeros(nx,nx);
A = eye(nx);

nf  = length(fStencil);
nfS = length(nfStencil);

Ad = zeros(nx, length(fStencil));
Bd = zeros(nx,length(nfStencil));

for i = 1:nx
    % Inserci�n del stencil en la malla para el nodo de estudio
    [fStencilMesh,nfStencilMesh] = locateStencil(i,nx,fStencil,nfStencil);
    % C�lculo de las distancias de la expansi�n de Taylor
    [fDistance,nfDistance0,s1,s2] = nodalDistance(fStencilMesh,nfStencilMesh,x,i);
    % Pivotaje respecto a f^{(n)}_{i} -> \alpha_i = 1. Por lo tanto se
    % elimina la informaci�n respectiva a dicho c�lculo
    nfDistance = nfStencilPivot(nfDistance0,i);
    % Taylor expansion
    [fMat,nfMat,lf,lnf] = taylorExpansion(fDistance,nfDistance,derivada);
    % Build LES matrix and check matrix condition
    [C,d] = buildMatrix(fMat,nfMat,derivada);
    % LES resolution
    coef = C\d;
    % Coef placement - Matrix A
    coefA = coef(1:lnf);
    for j1 = 1:lnf
        A(i,nfDistance(j1,2)) = coefA(j1);
    end
    % Coef placement - Matrix B
    coefB = coef(lnf+1:lf+lnf);
    for j2 = 1:lf
        B(i,fDistance(j2,2)) = coefB(j2);
    end
end

if flag==2 % extract diagonals
    for i=1:nf
        ind = i-(nf+1)/2;
        if ind<0 
            Ad(:,i) = [zeros(abs(ind),1);diag(A,ind)];
        else
            Ad(:,i) = [diag(A,ind);zeros(abs(ind),1)];
        end
    end
       
    for i=1:nfS
        ind = i-(nfS+1)/2;
        if ind<0 
            Bd(:,i) = [zeros(abs(ind),1);diag(B,ind)];
        else
            Bd(:,i) = [diag(B,ind);zeros(abs(ind),1)];
        end
    end 
end
        
    
    