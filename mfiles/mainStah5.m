%
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
% Script mainStah5. 
% 
% This routine reads all the statiscs and turbulent budgets of a pipe
% simulation
%
% Version 8.02.23

% Get the name of the host. 

[ret, name] = system('hostname');

if isunix
    root='/home/sergio/CFD-TFM Dropbox/pipe10k/';
    dataSpe = '/data/10k/spectraFull/';

        addpath(genpath(strcat(root,'mfiles/balh5'))); % h5 adapted budgets diles
        addpath(genpath(strcat(root,'mfiles/sta')));   % h5 adapted statistics diles
        addpath(genpath(strcat(root,'mfiles/general')));
        addpath(genpath(strcat(root,'mfiles/data')));
        addpath(genpath(strcat(root,'mfiles/figures')));
        addpath(genpath(strcat(root,'mfiles/derivatives')));
else
    if strcmp(name(1:min(5,length(name))),'SHEQF')
        root  = 'C:\Users\serhocal\CFD-TFM Dropbox\pipe10k\';

        addpath(genpath(strcat(root,'mfiles\balh5'))); % h5 adapted budgets diles
        addpath(genpath(strcat(root,'mfiles\sta')));   % h5 adapted statistics diles
        addpath(genpath(strcat(root,'mfiles\general')));
        addpath(genpath(strcat(root,'mfiles\data')));
        addpath(genpath(strcat(root,'mfiles\figures')));
        addpath(genpath(strcat(root,'mfiles\derivatives')));
    elseif strcmp(name(1:min(8,length(name))),'MarcosPM')
        root  = 'C:\Users\mpied\CFD-TFM Dropbox\pipe10k\';

        addpath(genpath(strcat(root,'mfiles\balh5'))); % h5 adapted budgets diles
        addpath(genpath(strcat(root,'mfiles\sta')));   % h5 adapted statistics diles
        addpath(genpath(strcat(root,'mfiles\general')));
        addpath(genpath(strcat(root,'mfiles\data')));
        addpath(genpath(strcat(root,'mfiles\figures')));
        addpath(genpath(strcat(root,'mfiles\derivatives')));

    else
        root = 'D:\sergio\CFD-TFM Dropbox\pipe10k\';        
    end
end






rootbal  = strcat(root,'/tbb/');
rootb    = strcat(root,'/binaries/');
rootsta  = strcat(root,'/sth/');
database = strcat(root,'/data/');
plotroot = strcat(root,'/figures/');

% Data fro plott


% Boundaries of the outer region
plotea    = 0; % if plotea = 1 we save figures to disk 


colors    = {'b','r','k','c-','b-','r:','g--','b*-','g:','b--'};
markers   = {'b-','r--','k-','c-','b-','r:','g--','b*-','g:','b--'};
markers2  = {'o','s','d','v','p','h','^'};

% Read data from 2000, 5000 and 4200. stored it in a variable.

%read180; % This file is in data. Reads the 180 from KTH

% caso = 'P125_21pir2';

datasets = {'pipe180lr'};
leyenda = datasets;
nf = length(datasets);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WHAT DO YOU want to do? 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Do we need to recompute statistis or turbulent budgets?
computeSta = 1; % 0 if you don't want to compute new stats
computeBal = 2;% 0 if you don't want to compute new stats, 2 if you don want to do anything at all

% Validation

validateBal = 0;
validateSta = true;

% errors

errorBars   = 0;
meaningBars = 0;

% moments
plotMomentsCenter = false;
plotMomentsLog   = false;
% balances
plotbal = false;
plotStat = false;


% Adding statistics

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Computing zone
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if computeSta==1
    for ii=1:nf
    sta2dat % Adds and adimensionalise the database. 
    sta(ii).stas=statTMP;
    sta(ii).case=datasets{ii};
    end
else
    for ii=1:nf
        load(strcat(rootb,datasets{ii},'_stat.sth.mat'),'statTMP');
        sta(ii).stas = statTMP;
    end
end

if computeBal==1
    frout = addbalh5(rootbal,datasets{1},rootb); % Suma la base de datos si frout=[];
    [stats,tbe]=varfitbal(frout);
    bal.stats = stats;
    bal.tbe = tbe;
    save(strcat(rootb,datasets{1},'_tbe.mat'),'bal');
elseif computeBal==2
    fprintf('Doing nothing...\n')
else
    load(strcat(rootb,datasets{1},'_tbe.mat'),'bal');
    signo = [1 -1 1 1 1 1];
    my2 = 1051;
    balvav = zeros(statTMP.my,6,6);
    balvavs = zeros(my2,6,6);
    I1 = 1:my2;
    I2 = statTMP.my:-1:my2;
    
    for ii=1:6
        balvav(:,1,ii) = bal.tbe.visd(:,ii);
        balvav(:,2,ii) = bal.tbe.prod(:,ii);
        balvav(:,3,ii) = bal.tbe.pdif(:,ii);
        balvav(:,4,ii) = bal.tbe.pstr(:,ii);
        balvav(:,5,ii) = bal.tbe.turb(:,ii);
        balvav(:,6,ii) = bal.tbe.disp(:,ii);
        
        balvavs(:,1,ii) = 0.5*(bal.tbe.visd(I1,ii) + signo(ii)*bal.tbe.visd(I2,ii));
        balvavs(:,2,ii) = 0.5*(bal.tbe.prod(I1,ii) + signo(ii)*bal.tbe.prod(I2,ii));
        balvavs(:,3,ii) = 0.5*(bal.tbe.pdif(I1,ii) + signo(ii)*bal.tbe.pdif(I2,ii));
        balvavs(:,4,ii) = 0.5*(bal.tbe.pstr(I1,ii) + signo(ii)*bal.tbe.pstr(I2,ii));
        balvavs(:,5,ii) = 0.5*(bal.tbe.turb(I1,ii) + signo(ii)*bal.tbe.turb(I2,ii));
        balvavs(:,6,ii) = 0.5*(bal.tbe.disp(I1,ii) + signo(ii)*bal.tbe.disp(I2,ii));
    end
    
end

% validation part

if validateBal
    validationBal;
end

if validateSta
    validationSta
end

if errorBars % Esto es solo para el error que usamos en la revisión.
    [y,vv,mom,nacum,head] = staOldWay(strcat(rootsta,datasets{1}));
    [yp2,tfrac,vav,vsq,mv,mvr,msq,msqr] = varfitstanew(y,vv,mom,nacum,head);
end

if plotMomentsCenter
    % Center
    momUCenter
    loglog_vw_centre
end
if plotMomentsLog

    % log layer
    momentLog
end


if plotbal
    fig1 = 1; % Viscous and buffer layer. Untill 200^+
    fig2 = 1; % Outer region: Whole region but premultiplied.
    fig3 = 0; % Plot pressure.
    escribe = 0; % write database
    balplot
end

if plotStat
    %figure1
    %figure2
    figure3
end
    





