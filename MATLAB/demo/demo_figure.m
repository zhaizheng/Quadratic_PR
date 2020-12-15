theta = pi*(rand([50,1])-0.5)*2;
x = sin(theta)+0.09*randn(size(theta));
y = cos(theta)+0.09*randn(size(theta));
plot(x,y,'*');
z = sqrt(x.^2+y.^2);
hold on
plot(x./z,y./z,'ro');
hold on
theta2 = 0:0.1:2*pi;
x1 = cos(theta2);
y1 = sin(theta2);
hold on
plot(x1,y1,'b:');

hold on

%h1 = quiver(x, y, (x./z-x), (y./z-y));
%set(h1,'AutoScale','off', 'AutoScaleFactor', 2)
%hold on


axis([-1.2 1.2 -1.2 1.2]);