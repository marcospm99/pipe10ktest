function [val,i,j] = min2D(A)

[VAL,JJ] = min(A);
[val,j] = min(VAL);
i=JJ(j);
