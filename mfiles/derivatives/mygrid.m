function [malla]=mygrid(N)

jj=0:N;   
dyeta=1.6; a=(.6*dyeta)^(4/3);
y1=.3; 

b=50; L=40;

c=((1+b).^(4/3)-b.^(4/3)-y1/a)./log((1+b)./b);
opts = odeset('reltol',1e-12);
[j,y]= ode45(@dydj,jj,0,opts,a,b,c,N,L);
dy=dydj(j,y,a,b,c,N,L);

retau=y(end), h=1/N;
y=y/y(end)-1;   dy=retau*h./dy;
y=[y; -y(end-1:-1:1)];  
dy=[dy; dy(end-1:-1:1)];
% save
% 
% a = 2*N+1;

% malla en asccii
malla = [y dy];

% save ('malla2.dat','malla','-ascii','-double')
% 
% dd=diff(y); dd=dd(1:N)*retau; 
% ym=(y(1:end-1)+y(2:end))/2;  ym=(1+ym(1:N))*retau;
% plot(ym,dd./0.8./ym.^(1/4))


% --------------------------------------
function dy=dydj(j,y,a,b,c,N,L)
N = N*1d0;

dyc=a*(4/3*(N+b).^(1/3)-c./(N+b));

dy=a*(4/3*(j+b).^(1/3)-c./(j+b));
dy=dyc+ (dy-dyc).*(1-1./cosh((j-N)/L));



