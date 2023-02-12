

dirinp = '/data/pipe/post/';
A=dir(strcat(dirinp,'*h5'));
for ii = 1:length(A)
    fprintf("Processing %i of %i \n",ii,length(A))
    frt = strcat(dirinp,A(ii).name);
    field2sta(frt)
end
