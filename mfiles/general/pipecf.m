function varargout=pipecf(A,i,varargin)
%
%--------------------------------------------------------------------------
% #         ###    #####  #######          #####    ###   #     # #     # -
% #          #    #     # #     #         #     #  #   #   #   #   #   #  -
% #          #    #       #     #               # # #   #   # #     # #   -
% #          #     #####  #     #          #####  #  #  #    #       #    -
% #          #          # #     #         #       #   # #   # #     # #   -
% #          #    #     # #     #         #        #   #   #   #   #   #  -
% #######   ###    #####  #######         #######   ###   #     # #     # -
%--------------------------------------------------------------------------
%
% function coucf
% This function plots the cf files of a LISO2000 simulation. To gain some
% flexibility it admits either a char A cointining the cf files to plot or
% a numerical A.
%
% Inputs
% A: Matrix containing the data to plot or the cf files. It could be empty.
%   In this case, it will plot all the cf files in the folder.
% i: input/output argument. if i is a scalar less than the size of A, the
% column i of A is plotted, and the output is i+1.
% On the contrary, if i is a scalar greater than 50 a group of plot are
% displayed. See below for details.
% if A is a character, i contains a struct array with all the data read.
% if i is a vector, the group of coulms of A given by i are plotted.
% varargin: more data, similar to A.
%
%

% Structure of the cf files.
labels={'\step'               ,... % 1
    '\Nr'            ,... % procs NR
    '\Ns'           ,... % procs Ns
    'time'   ,... % 4
    'timr'   ,... % 5
    'Ub'               ,... % 6
    'Uc'                 ,... % 7
    'Ufr'                 ,... % 8
    'Re_{\tau}'                 ,... % 9
    
};

if nargin==0
    A=[];
    i=-1;
end
ctrl = 0;
switch i
    case 93 % Thermal & momentum
        i = [6 4 5 19 21 22 3 14 16];
        lx=3;lz=3;
    case 99
        i = [6 7 8 9];
        ctrl = 1;
        lx=2;lz=2;
end


if nargin==1; varargout{1}=2; end

if (isempty(A)) || (ischar(A))
    if isempty(A)
        file=dir('*.cf');
    else
        file=dir(A);
    end
    
    if isempty(file)
        error('No cf file, exiting');
    end
    
    [n,~]=size(file);
    C=cell(n,2);
    for j =1:n
        A=load(file(j).name);
        C{j,1}=A;
        C{j,2}=file(j).name;
        %plot(A(:,1),A(:,i),strcat(lisocolors(j),lisosymbol(j)));
    end
    varargout{1}=C;
    
    if i<0
        return
    end
    
    figure(1); clf; hold on;
    if length(i)==1;
        for j =1:n;
            plot(C{j,1}(:,4),C{j,1}(:,i),strcat(lisocolors(j),lisosymbol(j)));
        end
        xlabel('t, code units')
        ylabel(labels{i});
    else
        for k=1:lx*lz
            subplot(lx,lz,k); hold on;
                for j =1:n;
                    plot(C{j,1}(:,4),C{j,1}(:,i(k)),strcat(lisocolors(j),lisosymbol(j)));
                end
            xlabel('t, code units')
            ylabel(labels{i(k)});
        end
    end
    
    if ctrl==3
        for j =1:n
            fprintf('file: %s, Average transposes %f, communications %f, time %f\n',C{j,2},mean(C{j,1}(:,18)),mean(C{j,1}(:,17)),mean(C{j,1}(:,16)))
        end
    end
    legend(C{:,2})
    return
end


%dt=mean(B(:,16));
%ds=mean(diff(B(:,1)))/5;
%vel=ds/dt;
%fprintf('Average time %i, Average Space %i, Average velocity %i \n',dt,ds,vel);
%     figure(2); clf; hold on;
%     re=2500;
%     v0=0.001;
%     uwall=1.789;
%     for j=1:n
%         dt=abs((C{j,1}(:,2)+C{j,1}(:,3)))/re;
%         plot(C{j,1}(:,1),dt,strcat(colors(mod(j,5)+1),symbols(mod(j,4)+1)));
%         plot(C{j,1}(:,1),v0*uwall)
%     end


[~,m] = size(A);

if (nargin>1) && (length(i)==1)
    if i>m; i=2;end
end

if nargin==1
    subplot(3,3,1); plot(A(:,1),A(:,3),'.')
    ylabel ('Utau')
    subplot(3,3,2); plot(A(:,1),A(:,4),'.')
    ylabel ('RE_{\tau}')
    subplot(3,3,3); plot(A(:,1),A(:,5),'.')
    ylabel ('ener(u)')
    subplot(3,3,4); plot(A(:,1),A(:,12),'.')
    ylabel ('ener(p)')
    subplot(3,3,5); plot(A(:,1),A(:,17),'.')
    ylabel ('Tiempo')
    subplot(3,3,6); plot(A(:,1),A(:,16),'.')
    ylabel ('Comunicaciones')
    subplot(3,3,7); plot(A(:,1),A(:,13),'.')
    ylabel ('mass(U)')
    subplot(3,3,8); plot(A(:,1),A(:,14),'.')
    ylabel ('mass(V)')
    subplot(3,3,9); plot(A(:,1),A(:,15),'.')
    ylabel ('Mass(W)')
    
    
elseif nargin==3
    B = varargin{1};
    plot(A(:,1),A(:,i),'.',B(:,1),B(:,i),'o')
elseif nargin==4
    B = varargin{1};
    C = varargin{2};
    plot(A(:,1),A(:,i),'.',B(:,1),B(:,i),'o',C(:,1),C(:,i),'*')
end

varargout{1}=i+1;
