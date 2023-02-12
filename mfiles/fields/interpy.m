function interpy(filinp,filout)
cd /data/pipe/fields
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
% This file interpolates the file in Fourier Space. The code reads M and
% compare it to the one of the field. If it is larger, the new file is
% paded with 0s.
% The order is fir K=6, M=4, for a total of (2*K-1)*m + k is 
%
% m\k | -5  -4  -3  -2  -1   0   1   2   3   4   5 
% ------------------------------------------------
%  0  |                      0   1   2   3   4   5
%  1  |  6   7   8   9  10  11  12  13  14  15  16
%  2  | 17  18  19  20  21  22  23  24  25  26  27
%  3  | 28  29  30  31  32  33  34  35  36  37  38 

% Adding modes in K means adding zeros at the beggining and end of this
% matrix 
% Adding modes in Z means adding zeros at the end of the matrix. 

% This file computes the statis
% Important. It doesn't change anything else. To change to a new file, use
% readh5 after using this one

% Reads data

if exist(filout, 'file') == 2
   delete(filout);
end


% Header
% The code will not take all this values into account. 

h5lt_makedataset(filout,'/Re'    ,1,0, 'double');
h5lt_makedataset(filout,'/alpha' ,1,0, 'double');
h5lt_makedataset(filout,'/time'  ,1,0,'double');

h5lt_makedataset(filout,'/dt'    ,1,0.01,'double');
h5lt_makedataset(filout,'/dtcfl' ,1,-1,'double');
h5lt_makedataset(filout,'/dtcor' ,1,-1,'double');

% New MEsh

r   = h5read(filinp,'/r');
Kin = h5read(filinp,'/K');
Min = h5read(filinp,'/M');
Nold = length(r);

i_N  = 192;
rr = newR(i_N);

N = length(rr);

h5lt_makedataset(filout,'/N' ,1,N,'int32');
h5lt_makedataset(filout,'/M' ,1,Min,'int32');
h5lt_makedataset(filout,'/K' ,1,Kin,'int32');
h5lt_makedataset(filout,'/Mp',1,1,'int32');
h5lt_makedataset(filout,'/r',double(length(rr)),rr,'double'); 

% Newfields 
fields = {'/Ur/Re','/Ur/Im','/Uz/Re','/Uz/Im','/Ut/Re','/Ut/Im'};
lv = (2*Kin-1)*(Min-1) + Kin;
tmpout = zeros(N,lv); 
for ii=1:6
    ext = fields{ii};
    h5create(filout,ext ,double([N,lv]),'Datatype','double');
    tmpin =  h5read(filinp,ext);
    for jj=1:lv
        tmpout(:,jj) = interp1(r,tmpin(:,jj),rr);
    end
    %if sum(isnan(tmpout)); 'hay nans', end
    h5write(filout,ext,tmpout)
    tmpout = 0*tmpout;
end
