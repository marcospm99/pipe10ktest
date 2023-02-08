function netcdf2h5(filinp,filout)
% script assign_old2h5
%
% Important. It doesn't change anything else. To change to a new file, use
% readh5 after using this one

if exist(filinp, 'file') ~= 2
   error([ 'The file ', filinp, ' does not exist']);
end
if exist(filout, 'file') == 2
   delete(filout);
end


time = ncreadatt(filinp,"/",'t');
Re   = ncreadatt(filinp,"/",'Re');
alp  = ncreadatt(filinp,"/",'alpha');

r     = ncread(filinp,"r") ;
dt    = ncread(filinp,"dt");
dtcor = ncread(filinp,"dtcor");
dtcfl = ncread(filinp,"dtcfl");

% Fields 

Ur = ncread(filinp,"Ur");
Uz = ncread(filinp,"Uz");
Ut = ncread(filinp,"Ut");

[~,n2,~] = size(Ur);

K = ncreadatt(filinp,"Ur","K");
M = ncreadatt(filinp,"Ur","M");
Mp= ncreadatt(filinp,"Ur","Mp");

h5create(filout,'/Re'  ,1,'Datatype','double');
h5create(filout,'/alp' ,1,'Datatype','double');
h5create(filout,'/time',1,'Datatype','double');

h5create(filout,'/dt'  ,1,'Datatype','double');
h5create(filout,'/dtcor' ,1,'Datatype','double');
h5create(filout,'/dtcfl',1,'Datatype','double');

h5write(filout,'/Re',Re);
h5write(filout,'/alp',alp);
h5write(filout,'/time',time);

h5write(filout,'/dt',dt);
h5write(filout,'/dtcor',dtcor);
h5write(filout,'/dtcfl',dtcfl);


nr = length(r);
h5create(filout,'/r',nr,'Datatype','double');
h5write(filout,'/r',r);

% Copy to new file.

h5create(filout,'/K',1,'Datatype','int32'); % Fourier en x. 
h5create(filout,'/M',1,'Datatype','int32'); % Fourier en theta
h5create(filout,'/N',1,'Datatype','int32'); % Puntos radiales
h5create(filout,'/Mp',1,'Datatype','int32'); % Puntos radiales

h5write(filout,'/K',int32(K));
h5write(filout,'/M',int32(M));
h5write(filout,'/N',int32(nr));
h5write(filout,'/Mp',int32(Mp)); % This sould be oalways one

h5create(filout,'/Ur/Re',double([nr,n2]),'Datatype','double');
h5create(filout,'/Ur/Im',double([nr,n2]),'Datatype','double');
h5create(filout,'/Uz/Re',double([nr,n2]),'Datatype','double');
h5create(filout,'/Uz/Im',double([nr,n2]),'Datatype','double');
h5create(filout,'/Ut/Re',double([nr,n2]),'Datatype','double');
h5create(filout,'/Ut/Im',double([nr,n2]),'Datatype','double');

h5write(filout,'/Ur/Re',Ur(:,:,1)); 
h5write(filout,'/Ur/Im',Ur(:,:,2));
h5write(filout,'/Uz/Re',Uz(:,:,1));
h5write(filout,'/Uz/Im',Uz(:,:,2));
h5write(filout,'/Ut/Re',Ut(:,:,1));
h5write(filout,'/Ut/Im',Ut(:,:,2));






