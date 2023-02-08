function stas = addSth(folder)
% This file computes the statis
% Important. It doesn't change anything else. To change to a new file, use
% readh5 after using this one


if nargin==0
    folder = 'C:\Users\serhocal\CFD-TFM Dropbox\pipe10k\sth\';
end

files = dir(strcat(folder,'*sth'));

stas.utau = 0.0338302;
stas.Retau = 179.3;

nf = length(files);
stas.r = h5read(strcat(files(1).folder,'/',files(1).name),'/r');
i_K = double(h5read(strcat(files(1).folder,'/',files(1).name),'/K'));
i_M = double(h5read(strcat(files(1).folder,'/',files(1).name),'/M'));
points = 1./(9*i_K*i_M);

nr = length(stas.r);

stas.URP = zeros(nr,nf);
stas.URZ = zeros(nr,nf);
stas.UTP = zeros(nr,nf);
stas.UZP = zeros(nr,nf);

stas.UR0 = zeros(nr,nf);
stas.UT0 = zeros(nr,nf);
stas.UZ0 = zeros(nr,nf);

stas.num = zeros(nf,1);
stas.time = zeros(nf,1);


for ii =1:nf
    
    frout = strcat(files(1).folder,'/',files(ii).name);
    
    stas.num(ii) = double(h5read(frout,'/num'));
    stas.time(ii)= double(h5read(frout,'/time'));
    
    stas.UR0(:,ii) = h5read(frout,'/mean_ur')./stas.num(ii)*points;
    stas.UT0(:,ii) = h5read(frout,'/mean_ut')./stas.num(ii)*points;
    stas.UZ0(:,ii) = h5read(frout,'/mean_uz')./stas.num(ii)*points;

    stas.URP(:,ii) = h5read(frout,'/stdv_ur')./stas.num(ii)*points;
    stas.UTP(:,ii) = h5read(frout,'/stdv_ut')./stas.num(ii)*points;
    stas.UZP(:,ii) = h5read(frout,'/stdv_uz')./stas.num(ii)*points;
    stas.URZ(:,ii) = h5read(frout,'/stdv_rz')./stas.num(ii)*points;
end

stas.urm = mean(stas.UR0,2);
stas.uzm = mean(stas.UZ0,2);
stas.utm = mean(stas.UT0,2);

uz0 = stas.uzm;%!-1d0+stas.r.^2;

stas.uzp = sqrt(mean(stas.UZP,2)-uz0.^2)./stas.utau;
stas.urz = mean(stas.URZ,2)/stas.utau.^2;
stas.urp = sqrt(mean(stas.URP,2))./stas.utau;
stas.utp = sqrt(mean(stas.UTP,2))./stas.utau;

stas.uzm = (stas.uzm+1d0-stas.r.^2)./stas.utau;

stas.r2 = (stas.r-1)*-1;

%tor180 = load('/home/sergio/CFD-TFM Dropbox/pipe10k/sth/Re180.prof');
tor180 =  load('C:\Users\serhocal\CFD-TFM Dropbox\pipe10k\sth\Re180.prof');
subplot(3,2,1); plot(stas.r2,stas.urz,tor180(:,1),-tor180(:,11),'o')
subplot(3,2,2); plot(stas.r2,stas.uzm,tor180(:,1),tor180(:,3),'o')
subplot(3,2,3); plot(stas.r2,stas.utm)
subplot(3,2,4); plot(stas.r2,stas.urp,tor180(:,1),tor180(:,5),'o')
subplot(3,2,5); plot(stas.r2,stas.uzp,tor180(:,1),tor180(:,4),'o')
subplot(3,2,6); plot(stas.r2,stas.utp,tor180(:,1),tor180(:,6),'o')




