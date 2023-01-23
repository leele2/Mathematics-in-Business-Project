function y = linpol(x0,y0,x1,y1,x)
    y = ( (y0*(x1-x)) + (y1*(x-x0)) ) / (x1-x0);
end