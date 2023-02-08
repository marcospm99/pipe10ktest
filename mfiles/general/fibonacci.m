function fb=fibonacci(my)
fb=zeros(my,1);
fb(1)=1;
fb(2)=1;
i=3;
while fb(i-1)<my
    fb(i) = fb(i-1)+fb(i-2);
    i=i+1;
end
fb([1 i-1:end])=[];