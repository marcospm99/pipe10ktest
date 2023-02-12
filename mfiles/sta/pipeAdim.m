function stats = pipeAdim(filin)
%--------------------------------------------------------------------------
% #         ###    #####  #######          #####    ###   #     # #     # -
% #          #    #     # #     #         #     #  #   #   #   #   #   #  -
% #          #    #       #     #               # # #   #   # #     # #   -
% #          #     #####  #     #          #####  #  #  #    #       #    -
% #          #          # #     #         #       #   # #   # #     # #   -
% #          #    #     # #     #         #        #   #   #   #   #   #  -
% #######   ###    #####  #######         #######   ###   #     # #     # -
%--------------------------------------------------------------------------
% Read .sth files
% 
stats.name = filin;

stats.nacum = double(h5read(filin,'/header/num'));

i_K= double(h5read(filin,'/header/K')); 
i_M= double(h5read(filin,'/header/M')); 

points = 1./(9*i_K*i_M);

stats.r = h5read(filin,'/header/r') ;
stats.timev = h5read(filin,'/header/timev');
stats.alpha = h5read(filin,'/header/alpha');

stats.urm = h5read(filin,'/sta/mean_ur')./stats.nacum*points;
stats.utm = h5read(filin,'/sta/mean_ut')./stats.nacum*points;
stats.uzm = h5read(filin,'/sta/mean_uz')./stats.nacum*points;

% stats.w1m = h5read(filin,'/sta/w1m');
% stats.w2m = h5read(filin,'/sta/w2m');
% stats.w3m = h5read(filin,'/sta/w3m');

% dU\Dy is w3m. 
stats.Re  = h5read(filin,'/header/Re');
stats.N   = h5read(filin,'/header/N');


stats.u_tau  = h5read(filin,'/sta/utau')./stats.nacum;
stats.u_cl  = h5read(filin,'/sta/ucl')./stats.nacum;
stats.Re_tau = stats.u_tau*stats.Re;


% Adimensionalise

fu = 1./stats.u_tau;
% fw = 1./(stats.u_tau^2*stats.Re)./stats.nacum;


stats.urp = sqrt(h5read(filin,'/sta/stdv_ur')./stats.nacum*points-stats.urm.^2)*fu;
stats.utp = sqrt(h5read(filin,'/sta/stdv_ut')./stats.nacum*points-stats.utm.^2)*fu;
stats.uzp = sqrt(h5read(filin,'/sta/stdv_uz')./stats.nacum*points-stats.uzm.^2)*fu;


% vorticies

%stats.w1p = sqrt(h5read(filin,'/sta/w1p')-stats.w1m.^2)*fw;
%stats.w2p = sqrt(h5read(filin,'/sta/w2p')-stats.w2m.^2)*fw;
%stats.w3p = sqrt(h5read(filin,'/sta/w3p')-stats.w3m.^2)*fw;

%stats.w1m = stats.w1m *fw;
%stats.w2m = stats.w2m *fw;
%stats.w3m = stats.w3m *fw;

% Reynolds stresses

stats.uruz  = (h5read(filin,'/sta/stdv_rz')./stats.nacum*points-stats.urm.*stats.uzm)*fu*fu;
%stats.uw = (h5read(filin,'/sta/uwr')-stats.um.*stats.wm)*fu*fu;
%stats.vw = (h5read(filin,'/sta/vwr')-stats.vm.*stats.wm)*fu*fu;



stats.urm = stats.urm*fu;
stats.utm = stats.utm*fu;
stats.uzm = (stats.uzm+1d0-stats.r.^2)*fu;

% Save y.



