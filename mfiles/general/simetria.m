function [ varargout ] = simetria( varargin )
%simetria Devuelve la mitad de un vector sim�trico haciendo la media
% En caso de ser antisim�trico tambi�n funciona, devolviendo el valor
% positivo.

% varargout = varargin;

for i=1:nargin
    vect = varargin{i};
    l = length(vect);
    j = floor(l/2);
    if mod(l,2) % Impar
        v1 = vect(j+1:-1:1);
        v2 = vect(j+1:end);
    else % Par
        v1 = vect(j:-1:1);
        v2 = vect(j+1:end);
    end
    
    if sum(v1.*v2) < 0 %Antisim�trico
        v = abs( v1-v2 )/2;
    else
        v =    ( v1+v2 )/2;
    end
    
    varargout{i} = v(end:-1:1);
    
end

