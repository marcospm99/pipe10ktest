function [frout,nf] = addsthPipe(rootin,filin,rootout)
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
root = strcat(rootin,filin,'/');

frout = [rootout,filin,'_avg.sth'];

if exist(frout,'file')~=0
    delete(frout)
end

files = dir(strcat(root,'*sth'));
nf = length(files);

% There are two groups. First we read headers and check that all of them
% corresponds to the same field.

ll=h5info(strcat(root,files(1).name));

header = ll.Groups(1).Datasets;

for i=1:length(header)
    
    if strcmp(header(i).Datatype.Class,'H5T_INTEGER')
        tmp='int32';
    elseif strcmp(header(i).Datatype.Class,'H5T_FLOAT')
        tmp='double';
    else
        error('check type!!!!')
    end
    
    h5create(frout,['/header/',header(i).Name],header(i).Dataspace.Size,'Datatype',tmp);
    
    ii = 1; 
    file = strcat(root,files(ii).name);
    
    tmp = h5read(file,['/header/',header(i).Name]);
    h5write(frout,['/header/',header(i).Name],tmp);
end

% Special cases: num, qav01, qav0L,time. 

num = 0;

timev = zeros(nf,1);

for ii = 1:nf
    file = strcat(root,files(ii).name);
    num = num + h5read(file,'/header/num');
    timev(ii) =  h5read(file,'/header/time');
end

h5write(frout,'/header/num/',num);
h5create(frout,'/header/timev',nf,'Datatype','double');
h5write(frout,'/header/timev/',timev);



% Statistics. 

sta =  ll.Groups(2).Datasets;

for i=1:length(sta)
    if strcmp(sta(i).Datatype.Class,'H5T_INTEGER')
        tmp='int32';
    elseif strcmp(sta(i).Datatype.Class,'H5T_FLOAT')
        tmp='double';
    else
        error('check type!!!!')
    end
    h5create(frout,['/sta/',sta(i).Name],sta(i).Dataspace.Size,'Datatype',tmp);
    
    if length(sta(i).Dataspace.Size) == 1
        tmp = zeros(sta(i).Dataspace.Size,1);
    else
        tmp = zeros(sta(i).Dataspace.Size);
    end
    
    for ii = 1:nf
        file = strcat(root,files(ii).name);
         try
        tmp = tmp + h5read(file,['/sta/',sta(i).Name]);
        catch 
             fprintf('problema con objeto %s en %s \n',sta(i).Name,file)
         end
    end
    
    % Save the value, not the mean

    h5write(frout,['/sta/',sta(i).Name],tmp);
end
