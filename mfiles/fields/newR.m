function r = newR(i_N)

N_ = i_N+floor((sqrt(i_N)));

for n = N_-i_N+1: N_
    r(n-N_+i_N,1) = 0.5d0*( 1d0+cos(pi*(N_-n)/N_) );
end
for n = 1: 10
    dr = 1.5d0*r(1)-0.5d0*r(2);
    r = r*(1d0+dr) - dr;
end