%--------------------------------------------------------------------------
%                                            #      ###                   -
%  #####      #    #####   ######           ##     #   #   #    #         -
%  #    #     #    #    #  #               # #    # #   #  #   #          -
%  #    #     #    #    #  #####             #    #  #  #  ####           -
%  #####      #    #####   #                 #    #   # #  #  #           -
%  #          #    #       #                 #     #   #   #   #          -
%  #          #    #       ######          #####    ###    #    #         -
%                                                                         -
%--------------------------------------------------------------------------
%
% function staMainthe
%
% This functions reads the raw data coming from one or several datasets. 
% For every dataset,
%   1. The function addsth averages each dataset. An _avg file is created
%   in h5 format. 
%   2. The function pipeSthreadh5 converts the raw data into useful one. A struct
%   array is created containing all the information needed
%   3. This array is save in a -mat format with a _stat extension
%
%   Last modification 10/02/23


if ~exist('datasets','var')
    error('database does not exist');
end


 for ind = 1:nf  % Read all datasets
     [frout,~] = addsthPipe(rootsta,datasets{ind},rootb); 
     statTMP = pipeAdim(frout);
     save(strcat(rootb,datasets{ind},'_stat.sth.mat'),'statTMP');
 end

