function fieldout = channelCenter(filin)
%
% Esta función hace la traspuesta local de datos escritos en la form
% ny,mz,mx a mx,mz,my
%
% Es cara, ojito a lso if para controlarla

S = fieldnames(filin);
ns = length(S);

for ii=1:ns
    tmp = filin.(S{ii});
    if S{ii}(1)~='j'
        tmp(:,:,end) = 2*tmp(:,:,end);
    end
    fieldout.(S{ii}) = tmp;
end

