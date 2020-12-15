x = -1:0.1:1;
x1 = [-1:0.2:-0.2,0.2:0.2:1];
x11 = -1:0.2:1;
y = x.^2;
y1 = x1.^2;
%plot(x,y,'linewidth',2);
%hold on
plot(x1,y1,'o','MarkerSize',8,'MarkerFaceColor','r');
hold on
z1 = x1.^2+0.4;
z = x.^2+0.4;
plot(x+0.1,z,'linewidth',2);
hold on
plot(x1+0.1,z1,'d','MarkerSize',8,'MarkerFaceColor','b');

a = 0:0.02:0.18;
b = 12.8*a.^2;
plot(a,b,'linewidth',1.5)
axis([-1 1 -0.2 1.4])
hold on
plot(0.18,12.8*0.18.^2,'d','MarkerSize',8,'MarkerFaceColor','m');
hold on
plot(0,0,'o','MarkerSize',8,'MarkerFaceColor','k');
hold on
plot(0.1,0.4,'d','MarkerSize',8,'MarkerFaceColor','k');
legend({'$\hat{\cal G}$','$\cal M$', '${\cal M}_{\hat {\cal G}}$','$\gamma(t)$','$\tilde{y}$','$\hat{x}_k$','$\tilde{y}_s$'},'Interpreter','latex')
set(gca,'FontSize',20);