function field=readFields(frt)
% script assign_old2h5
%
% Important. It doesn't change anything else. To change to a new file, use
% readh5 after using this one


if exist(frt, 'file') ~= 2
   error([ 'The file ', frt, ' does not exist']);
end

% Header
N = 384

for i=1:min(mx/2,mxe/2) % Number of planes.
    
    fprintf('Writting %i of %i \n',i,min(mx/2,mxe/2))
    fprintf('Leyendo\n')
    B = h5read(frt,'/vor',[1 1 double(i)],[double(2*my) double(mz) 1]);    
    np = single(interplane(B,my,mz,mze));
    h5write(frout,'/vor',np,[1 1 double(i)],[double(2*my) double(mze) 1]);
    B = h5read(frt,'/phi',[1 1 double(i)],[double(2*my) double(mz) 1]);
    np = single(interplane(B,my,mz,mze));
    h5write(frout,'/phi',np,[1 1 double(i)],[double(2*my) double(mze) 1]);
end


np = np*0;
for i = (min(mx/2,mxe/2) + 1):mxe/2
    h5write(frout,'/phi',np,[1 1 double(i)],[double(2*my) double(mze) 1]);
    h5write(frout,'/vor',np,[1 1 double(i)],[double(2*my) double(mze) 1]);
end


% save
%
function np=interplane(op,my,mz,mze)
np = zeros(2,my,mze);
op = reshape(op,2,my,mz);

mzm = min(mz,mze);

klen  = (mzm+1)/2;

np(:,:,1:klen)         = op(:,:,1:klen);
np(:,:,end-(0:klen-1)) = op(:,:,end-(0:klen-1));

np = reshape(np,2*my,mze);

function  np = interpy(B,y,ynew,my,mz)

%np = zeros(2,length(ynew),mz);
op = reshape(B,2,my,mz);
opt = permute(op,[2,1,3]);
np = interp1(y,opt,ynew);
np = reshape(permute(np,[2,1,3]),2*length(ynew),mz);




%
%
% function s=turdir(i,n)
%
% if i<0; display('negativo');s=[];return;end
% s=num2str(i);
%
% while length(s)<n
%     s=['0',s];
% end
%
% function [malla]=mygrid(N)
%
% jj=[0:N];
% dyeta=1.6; a=(.6*dyeta)^(4/3);
% y1=.3;
%
% b=50; L=40;
%
% c=((1+b).^(4/3)-b.^(4/3)-y1/a)./log((1+b)./b);
% opts = odeset('reltol',1e-12);
% [j,y]= ode45(@dydj,jj,0,opts,a,b,c,N,L);
% dy=dydj(j,y,a,b,c,N,L);
%
% retau=y(end), h=1/N;
% y=y/y(end)-1;   dy=retau*h./dy;
% y=[y; -y(end-1:-1:1)];
% dy=[dy; dy(end-1:-1:1)];
%
% malla = [y dy];
%
% % a = 2*N+1
% %
% % malla en asccii
% %
% % save ('malla2.dat','malla','-ascii','-double')
% %
% % dd=diff(y); dd=dd(1:N)*retau;
% % ym=(y(1:end-1)+y(2:end))/2;  ym=(1+ym(1:N))*retau;
% % plot(ym,dd./0.8./ym.^(1/4))
%
%
% % --------------------------------------
% function dy=dydj(j,y,a,b,c,N,L)
%
% dyc=a*(4/3*(N+b).^(1/3)-c./(N+b));
%
% dy=a*(4/3*(j+b).^(1/3)-c./(j+b));
% dy=dyc+ (dy-dyc).*(1-1./cosh((j-N)/L));
%
