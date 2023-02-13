%
%--------------------------------------------------------------------------
%                                            #      ###                    -
%  #####      #    #####   ######           ##     #   #   #    #         -
%  #    #     #    #    #  #               # #    # #   #  #   #          -
%  #    #     #    #    #  #####             #    #  #  #  ####           -
%  #####      #    #####   #                 #    #   # #  #  #           -
%  #          #    #       #                 #     #   #   #   #          -
%  #          #    #       ######          #####    ###    #    #         -
%                                                                         -
%--------------------------------------------------------------------------
%
% Script to validate pipes database. 

% Read data. comment and uncomment 
el180 = load(strcat(database,'180_Re_1.dat'));
leyenda = {'El Khoury et al','P180Low'};

figure; clf

% First plot. Plot comparison database

subplot(3,2,1); plot(el180(:,1),el180(:,3),'--g')
xlabel('1-r','fontsize',14);ylabel('$U_z^+$','fontsize',14,'interpreter','latex');

subplot(3,2,2); plot(el180(:,1),el180(:,5),'--g')
xlabel('1-r','fontsize',14);ylabel('$ur_{\rm rms}^+$','fontsize',14,'interpreter','latex');

subplot(3,2,3); plot(el180(:,1),el180(:,7),'--g')
xlabel('1-r','fontsize',14);ylabel('$uz_{\rm rms}^+$','fontsize',14,'interpreter','latex');

subplot(3,2,4); plot(el180(:,1),el180(:,6),'--g')
xlabel('1-r','fontsize',14);ylabel('$ut_{\rm rms}^+$','fontsize',14,'interpreter','latex');

subplot(3,2,5); plot(el180(:,1),el180(:,8),'--g')
xlabel('1-r','fontsize',14);ylabel('$uzur^+$','fontsize',14,'interpreter','latex');

plt = {'uzm','urp','uzp','utp','uruz'};

for ii=1:nf
    dt = (max(sta(ii).stas.timev)-min(sta(ii).stas.timev));
    utau = sta(ii).stas.u_tau;
    ucl = sta(ii).stas.u_cl;
    lt = 2./double(sta(ii).stas.alpha);
    fprintf('Sim %s. Time, eddy turnovers: %f, W.o.: %f \n', sta(ii).case, dt*utau,dt*ucl/(lt*pi))
    
    
    for jj=1:5
        r = (sta(ii).stas.r-1)*-1;
        tm = sta(ii).stas.(plt{jj});
        subplot(3,2,jj); hold on; plot(r,tm,markers{ii});
    end
end
    

subplot(3,2,1); legend(leyenda)

subplot(3,2,6); hold on

y = el180(:,1); 
[fStencil,nfStencil] = getStencil(5,5);
[Ax,Bx] = matrixCFD(y,fStencil,nfStencil,1,0);
dun = Ax\(Bx*el180(:,3));


ii=1;

r = (sta.stas.r-1)*-1;
[Ax,Bx] = matrixCFD(r,fStencil,nfStencil,1,0);
tm = sta.stas.uzm;
dtm = Ax\(Bx*tm);

plot(y*180,y.*dun,r*180,dtm.*r); set(gca,'xsca','lo')



% % Stencil de c�lculo
% [fStencil,nfStencil] = getStencil(5,5);
% %--------------------------------------------------------------------------
% %% Funci�n de c�lculo de las matrices A y B
% derivada = 1;
% [Ax,Bx] = matrixCFD(x,fStencil,nfStencil,derivada);
% 
% derivada = 2;
% [fStencil,nfStencil] = getStencil(7,7);
% [Axx,Bxx] = matrixCFD(x,fStencil,nfStencil,derivada);
% 
% %--------------------------------------------------------------------------
% f   =  sin(x);
% du  =  cos(x);
% ddu = -sin(x);
% dun = Ax\(Bx*f);
% ddun = Axx\(Bxx*f);







