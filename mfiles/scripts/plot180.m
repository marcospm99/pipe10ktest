% The three sets are assumed to be in the same folder

el180 = load('180_Re_1.dat');
p180l = load('pipe180.LR.dat');
p180d = load('pipe180.DR.dat');


subplot(3,2,1); plot(p180l(:,1),p180l(:,3),'k',p180d(:,1),p180d(:,3),'.r',el180(:,1),el180(:,3),'--g')
xlabel('1-r','fontsize',14);ylabel('$U_z^+$','fontsize',14,'interpreter','latex');
leyenda = {'P180Low','P180High','El Khoury et al'};
legend(leyenda);
subplot(3,2,2); plot(p180l(:,1),p180l(:,4),'k',p180d(:,1),p180d(:,4),'.r',el180(:,1),el180(:,5),'--g')
xlabel('1-r','fontsize',14);ylabel('$ur_{\rm rms}^+$','fontsize',14,'interpreter','latex');
subplot(3,2,3); plot(p180l(:,1),p180l(:,5),'k',p180d(:,1),p180d(:,5),'.r',el180(:,1),el180(:,7),'--g')
xlabel('1-r','fontsize',14);ylabel('$uz_{\rm rms}^+$','fontsize',14,'interpreter','latex');

subplot(3,2,4); plot(p180l(:,1),p180l(:,6),'k',p180d(:,1),p180d(:,6),'.r',el180(:,1),el180(:,6),'--g')
xlabel('1-r','fontsize',14);ylabel('$ut_{\rm rms}^+$','fontsize',14,'interpreter','latex');

subplot(3,2,5); plot(p180l(:,1),p180l(:,7),'k',p180d(:,1),p180d(:,7),'.r',el180(:,1),el180(:,8),'--g')
xlabel('1-r','fontsize',14);ylabel('$uzur^+$','fontsize',14,'interpreter','latex');




pl = diff(p180l(:,2));
pd = diff(p180d(:,2));
el = diff(el180(:,2));
subplot(3,2,6); plot(p180l(1:end-1,1),pl,'k',p180d(1:end-1,1),pd,'.r',el180(1:end-1,1),abs(el),'--g')
xlabel('1-r','fontsize',14);ylabel('$uzur^+$','fontsize',14,'interpreter','latex');
