clear all
clc
tic;
%% SCRIPT PRINCIPAL - ESQUEMA CFD 2D ACOPLADO
nx = 31;
ny = 31;
%-----------------
ax = 0;
bx = 1;
ay = 0;
by = 1;
%-----------------
x = getDiscretization(ax,bx,nx);
y = getDiscretization(ay,by,ny);
%% Malla de estudio
[yDfc,xDfc] = meshgrid(y,x);
%% Construcci�n de las matrices del sistema
[G,Ax,Ay] = matrixBuilderCFDTwoD(x,y);
%% Condiciones de contorno: DIRICHLET
uCCstep1 = xDfc.^2.*yDfc+yDfc.^2.*xDfc-3;
uCCOnes  = [ones(1,ny);ones(nx-2,1),zeros(nx-2,ny-2),ones(nx-2,1);ones(1,ny)];
uCC      = uCCstep1.*uCCOnes;
uCCv     = reshape(uCC,nx*ny,1);
%% Funci�n "f": u''=f
f = 2*(xDfc+yDfc);
%% Construcci�n de la matriz H: [G]�{u}={H}
H = reshape(Ax*f*Ay,nx*ny,1)-G*uCCv;
%% Incorporaci�n de las condiciones de contorno - influencia G y H
rowsNoNull = find(uCCOnes);
G(rowsNoNull,:) = [];
G(:,rowsNoNull) = [];
H(rowsNoNull,:) = [];
%% Condicionamiento G y H
[rowsG,~] = size(G);
GC = zeros(size(G));
HC = zeros(rowsG,1);
for i=1:rowsG
    condV1 = find(G(i,:));
    condV2 = sort(abs(G(i,condV1)));
    condH  = abs(H(i));
    cond   = min(condH,condV2(1));
    GC(i,:) = G(i,:)/cond;
    HC(i,:) = H(i,:)/cond;
end
%% Resoluci�n del SEL
u    = GC\HC;
uDFC = reshape(u,nx-2,ny-2);
u    = uCC+[zeros(1,ny);zeros(nx-2,1),uDFC,zeros(nx-2,1);zeros(1,ny)];
toc
%% Representaci�n gr�fica
plotResults = 1;
if plotResults
    %--
    figure(1)
    xlim([min(y),max(y)]);
    ylim([min(x),max(x)]);
    surf(yDfc,xDfc,u)
    set(gca,'FontSize',18);
    xlabel('Eje $Y$','Interpreter','Latex','FontSize',22);
    ylabel('Eje $X$','Interpreter','Latex','FontSize',22);
    zlabel('$u(x,y)$','Interpreter','Latex','FontSize',22);
    title('$u(x,y)=x^2y+xy^2-3$','Interpreter','Latex','FontSize',24)
    %--
    figure(2)
    surf(yDfc,xDfc,abs((u-uCCstep1)./uCCstep1))
    set(gca,'zScale','Log','FontSize',18);
    title('Error relativo - $\varepsilon_R$','Interpreter','Latex','FontSize',24)
    xlabel('y','Interpreter','Latex','FontSize',22);
    ylabel('x','Interpreter','Latex','FontSize',22);
    zlabel('$\varepsilon_R$','Interpreter','Latex','FontSize',22);
end