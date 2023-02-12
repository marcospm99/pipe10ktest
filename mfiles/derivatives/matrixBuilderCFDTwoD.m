function[G,Ax,Ay] = matrixBuilderCFDTwoD(x,y)
%--------------------------------------------------------------------------
% 8-MAYO-2014
%%=============================== IN
% x,y       -> vector con los puntos
% fpp_{x,y} -> número de puntos de segunda derivada (x,y)
% f_{x,y}   -> número de puntos de función
%%=============================== END IN
%%=============================== OUT
% [G]   -> Matriz correspondiente a la igualdad: B^x·u·A^y+A^x·u·B^y = G·u
% [A_x] -> Matriz de coeficientes de la segunda derivada corresp. a x
% [A_y] -> Matriz de coeficientes de la segunda derivada corresp. a y
%%=============================== END OUT
%--------------------------------------------------------------------------
%%============================= MATRICES
[Ax,Bx,AyNT,ByNT,fx,fppx,fy,fppy] = buildTwoDMatrices(x,y);%NT->No Transposed
%%============================= FIN MATRICES
%--------------------------------------------------------------------------
%%============================= MATRICES Y
Ay          = transpose(AyNT);
By          = transpose(ByNT);
%%============================= FIN MATRICES Y
%--------------------------------------------------------------------------
%%============================= MATRIZ G
p     = length(x);
q     = length(y);
nodes = p*q;
G1    = zeros(nodes);
G2    = zeros(nodes);
for k = 1:q
    for h = 1:p
        m1Limit = sumLimit(k,fppy,q);
        for m1 = m1Limit(1,1):m1Limit(1,2)
            a1Limit = sumLimit(h,fx,p);
            for a1 = a1Limit(1,1):a1Limit(1,2)
               G1(((k-1)*p+h),((m1-1)*p+a1)) = G1(((k-1)*p+h),((m1-1)*p+a1)) + Bx(h,a1)*Ay(m1,k);
            end
        end
        m2Limit = sumLimit(k,fy,q);
        for m2 = m2Limit(1,1):m2Limit(1,2)
            a2Limit = sumLimit(h,fppx,p);
            for a2 = a2Limit(1,1):a2Limit(1,2)
               G2(((k-1)*p+h),((m2-1)*p+a2)) = G2(((k-1)*p+h),((m2-1)*p+a2)) + Ax(h,a2)*By(m2,k);
            end
        end
    end
end
G=G1+G2;
%%============================= FIN MATRIZ G