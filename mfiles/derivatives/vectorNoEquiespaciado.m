function x = vectorNoEquiespaciado(puntos,p_end)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Esta funci�n devuelve un  vector de puntos no equiespaciado entre los
%puntos 0 y 2*p_end, estos inclu�dos. El n�mero de puntos debe ser impar. 
%Si no lo es, el programa lo har� impar
% Alfonso Pallar�s - Last check - 23-Sept-2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Comprobamos que es impar el n�mero de puntos
if mod(puntos,2)==0
    puntos = puntos + 1;
end
% %El espaciado se realiza con respecto a una ley exponencial y = e^x
% aux = linspace(0,log(1+p_end),puntos/2);
% x1  = exp(aux)-1;
% x2  = sort(2*p_end-x1);
% x   = transpose([x1,x2(1,2:end)]);

[malla]=mygrid((puntos-1)/2);
x=malla(:,1)

