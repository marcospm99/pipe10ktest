function interpxz(filinp,filout,M,K)
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
h5lt_makedataset(filout,'/time'  ,1,-1,'double');

h5lt_makedataset(filout,'/dt'    ,1,-1,'double');
h5lt_makedataset(filout,'/dtcfl' ,1,-1,'double');
h5lt_makedataset(filout,'/dtcor' ,1,-1,'double');



% New MEsh

r   = h5read(filinp,'/r');
Kin = h5read(filinp,'/K');
Min = h5read(filinp,'/M');
N = length(r);

h5lt_makedataset(filout,'/N' ,1,N,'int32');
h5lt_makedataset(filout,'/M' ,1,M,'int32');
h5lt_makedataset(filout,'/K' ,1,K,'int32');
h5lt_makedataset(filout,'/Mp',1,1,'int32');
h5lt_makedataset(filout,'/r',double(N),r,'double'); 

% Newfields 
fields = {'/Ur/Re','/Ur/Im','/Uz/Re','/Uz/Im','/Ut/Re','/Ut/Im'};
tmpout = zeros(N,(2*K-1)*M + K); 
for ii=1:6
    ext = fields{ii};
    h5create(filout,ext ,double([N, (2*K-1)*M + K]),'Datatype','double');
    tmpin =  h5read(filinp,ext);
    % first step 
    tmpout(:,1:Kin) = tmpin(:,1:Kin);
    ind1 = K;
    ind2 = Kin;
    vc = 1:2*Kin-1;
    
    for mm=1:Min-1
        pos1 = ind1 + K-Kin;% We start in the position of the mode m. There are 2*K-1 things here
        
        tmpout(:,pos1 + vc) = tmpin(:,ind2+vc);
        ind1 = ind1 + 2*K-1;
        ind2 = ind2 + 2*Kin-1;
    end
    h5write(filout,ext,tmpout)
    tmpout = 0*tmpout;
end
