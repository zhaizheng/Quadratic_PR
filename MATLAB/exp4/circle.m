function data = circle(n,sigma)
    %n = 1000 sigma = 0.2
    t = linspace(pi/2,5*pi,n);
    r = linspace(1,7,n);
    
    x = r.*cos(t)+sigma*randn(1,n);
    y = r.*sin(t)+sigma*randn(1,n);
    data = [x;y];
%     plot(x,y,'.')
end
