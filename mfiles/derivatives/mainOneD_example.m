clear all
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
n   = 101; % $n$ debe ser impar si se emplea la ley de discretizaci�n exp.
x   = getDiscretization(x_1,x_n,n,1);
nx  = length(x); % Por seguridad
%--------------------------------------------------------------------------
% Stencil de c�lculo
[fStencil,nfStencil] = getStencil(5,5);
%--------------------------------------------------------------------------
%% Funci�n de c�lculo de las matrices A y B
derivada = 1;
[Ax,Bx] = matrixCFD(x,fStencil,nfStencil,derivada);

derivada = 2;
[fStencil,nfStencil] = getStencil(7,7);
[Axx,Bxx] = matrixCFD(x,fStencil,nfStencil,derivada);

%--------------------------------------------------------------------------
f   =  sin(x);
du  =  cos(x);
ddu = -sin(x);
dun = Ax\(Bx*f);
ddun = Axx\(Bxx*f);
%--------------------------------------------------------------------------

% Resultados
errRel = abs((dun-du)./du);
%-----------------
figure(1)
plot(x,errRel,'LineWidth',2)
set(gca,'YScale','Log','FontSize',18)
title('Evoluci\''on del $\varepsilon_{R}$','Interpreter','LaTeX','FontSize',20)
xlabel('$\vec{x}$','Interpreter','LaTeX','FontSize',20)
ylabel('$\varepsilon_{R}$','Interpreter','LaTeX','FontSize',20)
%-----------------
figure(2)
plot(x,du,'o',x,dun,'*')
set(gca,'FontSize',18)
title('Resoluci\''on num\''erica','Interpreter','LaTeX','FontSize',20)
xlabel('$\vec{x}$','Interpreter','LaTeX','FontSize',20)
ylabel('$f(x)$','Interpreter','LaTeX','FontSize',20)
%-----------------
figure(3)
semilogy(x,abs((ddu-ddun)./ddun),'*')
set(gca,'FontSize',18)
title('Resoluci\''on num\''erica','Interpreter','LaTeX','FontSize',20)
xlabel('$\vec{x}$','Interpreter','LaTeX','FontSize',20)
ylabel('$f(x)$','Interpreter','LaTeX','FontSize',20)
