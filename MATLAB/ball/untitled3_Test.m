x = linspace(0.1, 10, 200);
a = 1;
y = x.^2.*exp(-(x./a).^2);
figure
plot(x,y,'-')