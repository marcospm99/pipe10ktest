function h5lt_makedataset(filout,str,dim,value,datatype)

if length(value)~=dim
    error('Check Dimensions')
end

h5create(filout,str ,dim,'Datatype',datatype);
h5write(filout,str ,value);