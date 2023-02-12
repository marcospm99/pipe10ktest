% clear all
%clc
%% MAIN SCRIPT - COMPACT FINITE DIFFERENCES - ONE DIMENSION
% Alfonso Pallar�s - 23-Sept-2014
%--------------------------------------------------------------------------
%% Preparaci�n del c�lculo
%--------------------------------------------------------------------------
% Construcci�n del vector discreto del dominio de estudio: [x_1,x_n]\in
% \mathbb{R}^{1} \rightarrow {x_1,x_n} \equiv {x_i}_{i=1}^{n}
x_1 = -1;
x_n = 1;
n   = 1001; % $n$ debe ser impar si se emplea la ley de discretizaci�n exp.
x   = getDiscretization(x_1,x_n,n,1);
nx  = length(x); % Por seguridad
%--------------------------------------------------------------------------
% Stencil de c�lculo
[fStencil,nfStencil] = getStencil(7,5);
%--------------------------------------------------------------------------
% Calculamos la derivada con una matriz equiespaciada

derivada = 1;
[Ax,Bx] = matrixCFD(x,fStencil,nfStencil,derivada);

%--------------------------------------------------------------------------

f    =  x.*sin(x);
due  =  sin(x)+x.*cos(x);
duen = Ax\(Bx*f);

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% Calculamos la derivada con una matriz no equiespaciada
xn   = getDiscretization(x_1,x_n,n,2);
[Axn,Bxn] = matrixCFD(xn,fStencil,nfStencil,derivada);

fn    =  xn.*sin(xn);
dune  =  sin(xn)+xn.*cos(xn);
dunen = Axn\(Bxn*fn);

fmap  = 1./(Ax\(Bx*xn));
malla = mygrid((n-1)/2);

fmap=malla(:,2);

duenmap = (Ax\(Bx*fn)).*fmap;


% Resultados
errRele  = abs((duen-due)./due);

errRelen = abs((dunen-dune)./dune);
errRelemap  = abs((duenmap-dune)./dune);

% errRelen = abs((duen-duen./fmap)./duen);

%-----------------
figure(1)
plot(x,errRele,'o',xn,errRelen,'*',xn,errRelemap,'d')
set(gca,'YScale','Log','FontSize',18)
title('Evoluci\''on del $\varepsilon_{R}$','Interpreter','LaTeX','FontSize',20)
xlabel('$\vec{x}$','Interpreter','LaTeX','FontSize',20)
ylabel('$\varepsilon_{R}$','Interpreter','LaTeX','FontSize',20)