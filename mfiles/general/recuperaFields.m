function fieldout = transposeFields(filin)
%
% Esta función hace la traspuesta local de datos escritos en la form
% ny,mz,mx a mx,mz,my
%
% Es cara, ojito a lso if para controlarla

S = fieldnames(filin);
ns = length(S);

for ii=1:ns
    fieldout.(S{ii}) = permute(filin.(S{ii}),[3,2,1]);
end

