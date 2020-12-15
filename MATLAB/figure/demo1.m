t = -1:0.2:1;
y = sqrt(1- 1/3*t.^2)+0.05*randn(size(t));
plot(t,sqrt(1- 1/3*t.^2),'--','LineWidth',1)
hold on
plot(t, y, '*','MarkerSize',8)
hold on
plot(0,0.9,'o', 'MarkerSize',8,'MarkerFaceColor','r')
t = 0:0.1:pi;
x = 0.8*cos(t);
y = 0.8*sin(t)+0.3;
plot(x,y,'-','LineWidth',2)
axis([-1, 1.1,0.1,1.3])