% function field2staidea2(frt)
% This file computes the statis
% Important. It doesn't change anything else. To change to a new file, use
% readh5 after using this one
% nh = (2*512-1)*127+512


frt = 'D:\calor\pipe180.1098.h5';

% Orden ur (modos k;modos M) 

frout = strcat(frt(1:end-2),'sta');

if exist(frt, 'file') ~= 2
   error([ 'The file ', frt, ' does not exist']);
end
if exist(frout,'file')~=0
    delete(frout)
end

% Header
time = h5read(frt,'/time');
Re   = h5read(frt,'/Re');
alp  = h5read(frt,'/alpha'); 

r = h5read(frt,'/r');

uz1 = h5read(frt,'/Ur/Re');
uz2 = h5read(frt,'/Ur/Im'); 

uz = uz1 + 1i * uz2; 

uz0 = uz(:,1);% +1d0-r.^2; % Ojito, hay que quitarle esto, por ahora un arcano

uzp = 2*sum(uz.*conj(uz),2) - sum(uz(:,1:512).*conj(uz(:,1:512)),2);

plot(sqrt(uzp)./0.0338)


