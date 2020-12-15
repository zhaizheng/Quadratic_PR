function base = base_Calulate(x,y)
    base = [x.^4, x.^3*y, x.^2*y.^2, x*y.^3, y.^4, x.^3, x.^2*y, x*y.^2, y.^3, x.^2, x*y, y.^2, x, y, 1];
end