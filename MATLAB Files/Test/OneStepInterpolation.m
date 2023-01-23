clear
%% Parameters
r = 0.3;
u  = 1.1;
d  = 1/u;
delta_t = 0.1;
p = (exp(r*delta_t) - d)/(u-d);
disf = exp(-r*delta_t);
%% Node info
c_node = 11.05;
u_node = 12.1;
d_node = 10;
Ku  = (c_node*2+u_node)/3;
Kd  = (c_node*2+d_node)/3;
avgs = [12.21, 11.05, 10, 9.11, 8.19];
opts = [2.21, 1.05, 0, 0, 0];
%% One-Step Calculation using intepolation
err = 1e-3;
[loc, ubound, lbound] = findInSorted(Ku,avgs,err);
if loc > 0
    Cu = opts(loc);
else
    Cu = linpol(avgs(lbound),opts(lbound),avgs(ubound),opts(ubound),Ku);
end
[loc, ubound, lbound] = findInSorted(Kd,avgs,err);
if loc > 0
    Cd = opts(loc);
else
    Cd = linpol(avgs(lbound),opts(lbound),avgs(ubound),opts(ubound),Kd);
end
out = disf*(p*Cu + (1-p)*Cd)