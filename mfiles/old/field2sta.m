function field2sta(frt)
% This file computes the statis
% Important. It doesn't change anything else. To change to a new file, use
% readh5 after using this one

frout = strcat(frt(1:end-2),'sta');

if exist(frt, 'file') ~= 2
   error([ 'The file ', frt, ' does not exist']);
end
if exist(frout,'file')~=0
    delete(frout)
end
r = h5read('/data/pipe/post/base/pipe180.1000.h5','/r');
% Header
time = h5read(frt,'/time');
Re   = h5read(frt,'/Re');
alp  = h5read(frt,'/alpha'); 

ur = h5read(frt,'/Ur/Re') + 1i * h5read(frt,'/Ur/Im'); 

uz = h5read(frt,'/Uz/Re') + 1i * h5read(frt,'/Uz/Im'); 

ut = h5read(frt,'/Ut/Re') + 1i * h5read(frt,'/Ut/Im'); 


uz0 = uz(:,1)+1d0-r.^2; % Ojito, hay que quitarle esto, por ahora un arcano
ur0 = ur(:,1);
ut0 = ut(:,1);

uzp = 2*sum(uz.*conj(uz),2) - sum(uz(:,1:512).*conj(uz(:,1:512)),2);
urp = 2*sum(ur.*conj(ur),2) - sum(ur(:,1:512).*conj(ur(:,1:512)),2);
utp = 2*sum(ut.*conj(ut),2) - sum(ut(:,1:512).*conj(ut(:,1:512)),2);

urz = real(2*sum(uz.*conj(ur),2) - sum(uz(:,1:512).*conj(ur(:,1:512)),2));

h5create(frout,'/Re'    ,1,'Datatype','double');
h5create(frout,'/alpha'   ,1,'Datatype','double');
h5create(frout,'/time'  ,1,'Datatype','double');

h5write(frout,'/Re',Re);
h5write(frout,'/time',time);
h5write(frout,'/alpha',alp);
% 00 modes

my = length(r);

h5create(frout,'/ur0' ,double(my),'Datatype','double');
h5create(frout,'/ut0' ,double(my),'Datatype','double');
h5create(frout,'/uz0' ,double(my),'Datatype','double');
h5create(frout,'/r'   ,double(my),'Datatype','double');

h5write(frout,'/ur0',ur0);
h5write(frout,'/uz0',uz0);
h5write(frout,'/ut0',ut0);

h5write(frout,'/r',r);

h5create(frout,'/urp' ,double(my),'Datatype','double');
h5create(frout,'/utp' ,double(my),'Datatype','double');
h5create(frout,'/uzp' ,double(my),'Datatype','double');
h5create(frout,'/urz' ,double(my),'Datatype','double');

h5write(frout,'/urp',urp);
h5write(frout,'/uzp',uzp);
h5write(frout,'/utp',utp);
h5write(frout,'/urz',urz);

movefile(frt,'/data/pipe/post/pasadas/')
movefile(frout,'/home/sergio/CFD-TFM Dropbox/pipe10k/sta/')

