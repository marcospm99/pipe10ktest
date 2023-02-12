function createDatabase(sta)
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
%
% This file creates an ascii file containing he statistics of the
% simulation
% All folders and data are in the sta dield
%
% 


FID = fopen([sta.data,'_DR.dat'],'w+');

fprintf(FID, '%s\n','%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
fprintf(FID, '%s\n','%---------------------------------------------------------------------%');
fprintf(FID, '%s\n','%                                       #      ###                   -%');
fprintf(FID, '%s\n','%  #####      #    #####   ######      ##     #   #   #    #         -%');
fprintf(FID, '%s\n','%  #    #     #    #    #  #          # #    # #   #  #   #          -%');
fprintf(FID, '%s\n','%  #    #     #    #    #  #####        #    #  #  #  ####           -%');
fprintf(FID, '%s\n','%  #####      #    #####   #            #    #   # #  #  #           -%');
fprintf(FID, '%s\n','%  #          #    #       #            #     #   #   #   #          -%');
fprintf(FID, '%s\n','%  #          #    #       ######     #####    ###    #    #         -%');
fprintf(FID, '%s\n','%                                                                    -%');
fprintf(FID, '%s\n','%---------------------------------------------------------------------%');
fprintf(FID, '%s\n','%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
fprintf(FID, '%s\n','%%%%%%                      2023 - 2024                          %%%%%%');
fprintf(FID, '%s\n','%%%%%%                                                           %%%%%%');
fprintf(FID, '%s\n','%%%%%%  Contact: sergio.hoyas@mot.upv.es                         %%%%%%');
fprintf(FID, '%s\n','%%%%%%                                                           %%%%%%');
fprintf(FID, '%s\n','%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
fprintf(FID, '\n');
fprintf(FID, '%s\n','%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
fprintf(FID, '%s\n','%%%%%%                                                           %%%%%%');
fprintf(FID, '%s\n','%%%%%%                     <<< CAUTION >>>                       %%%%%%');
fprintf(FID, '%s\n','%%%%%%                                                           %%%%%%');
fprintf(FID, '%s\n','%%%%%%  All rights are reserved by the authors.                  %%%%%%');
fprintf(FID, '%s\n','%%%%%%  No part of the data described herein may be represented  %%%%%%');
fprintf(FID, '%s\n','%%%%%%  without reference. The data base may be used without     %%%%%%');
fprintf(FID, '%s\n','%%%%%%  notification to the author''s laboratory, asuming that    %%%%%%');
fprintf(FID, '%s\n','%%%%%%  it has been properly cited.                              %%%%%%');
fprintf(FID, '%s\n','%%%%%%                                                           %%%%%%');
fprintf(FID, '%s\n','%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
fprintf(FID, '\n');
fprintf(FID, '%s\n','%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
fprintf(FID, '%s\n','%%                                                                   %%');
fprintf(FID, '%s\n','%%  NOMENCLATURE:                                                    %%');
fprintf(FID, '%s\n','%%   z = streamwise coordinate                                       %%');
fprintf(FID, '%s\n','%%   r = wall-normal coordinate                                      %%');
fprintf(FID, '%s\n','%%   t = Azimuthal coordinate                                        %%');
fprintf(FID, '%s\n','%%   R = Radius                                                      %%');
fprintf(FID, '%s\n','%%                                                                   %%');
fprintf(FID, '%s\n','%%                                                                   %%');
fprintf(FID, '%s\n','%%  BOX SIZE in R (z,r,t): 10pi|1|2pi                                %%');
fprintf(FID, '%s','%%  FLOW CONDITIONS:  Retau  = ');
fprintf(FID, '%3.0f ,',sta.Retau);
fprintf(FID, '%s\n%s','%%','%%');
fprintf(FID, '                     utau  = %f                              ',sta.utau);
fprintf(FID, '%s\n%s','%%','%%');
fprintf(FID, '%s\n','%%');
fprintf(FID, '%s\n','%%                                                                   %%');
fprintf(FID, '%s\n','%%    Data file created: FEB 2023                                    %%');
fprintf(FID, '%s\n','%%                                                                   %%');
fprintf(FID, '%s\n','%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
fprintf(FID, '\n');
fprintf(FID, '%%   (1-r)     ((1-r)^+)    (Uz^+) (Uz/U_b) (ur_{rms}^+) (ut_{rms}^+) (uz_{rms}^+) (uzur^+) ');
fprintf(FID,'\n');

m = length(sta.r);

for i=m:-1:1
    fprintf(FID,' %13.8f ', 1-sta.r(i));
    fprintf(FID,' %13.8f ', (1-sta.r(i))*sta.Retau);
    fprintf(FID,' %13.8f ', sta.uzm(i)); 
    fprintf(FID,' %13.8f ', sta.urp(i)); 
    fprintf(FID,' %13.8f ', sta.uzp(i)); 
    fprintf(FID,' %13.8f ', sta.utp(i)); 
    fprintf(FID,' %13.8f ', sta.urz(i));
    fprintf(FID,'\n');
end

fclose(FID);

%oldfile = [root_save,file,'.base'];
%newfile = [root_save,file,'.tkin'];

%copyfile(oldfile,newfile);
