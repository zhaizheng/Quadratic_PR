clear all;
format short;
syms x
u=[-1/4;-1/4;-1/4;3/4];
%u=[1;0;0;0];
t=4/3;
q1=i;q2=1+i;q3=1;q4=0;
q=[q1;q2;q3;q4];
k=5;
tic;
for n=1:k
    u1=[];
    t=min(double(vpasolve(36*x^3-48*x^2+15*x-t==0)))
    a=(9-42*t+36*t^2)/(3*(4-29*t+60*t^2-36*t^3));
    b=(6-24*t+18*t^2)/(3*(4-29*t+60*t^2-36*t^3));
    c=1/(3*(4-29*t+60*t^2-36*t^3)); 
    d=(2-3*t)/(3*(4-29*t+60*t^2-36*t^3));
    for j=1:4:length(u)
        u2=[1,0,0,0;
         a,c,c,c;
         b,d,d,d;
         a,c,c,c;
         c,a,c,c;
         0,1,0,0;
         c,a,c,c;
         d,b,d,d;
         d,d,b,d;
         c,c,a,c;
         0,0,1,0;
         c,c,a,c;
         c,c,c,a;
         d,d,d,b;
         c,c,c,a;
         0,0,0,1;
         b,d,d,d;
         d,b,d,d;
         d,d,b,d;
         d,d,d,b]*[u(j);u(j+1);u(j+2);u(j+3)];%+[0; (1-a-3*c)/4; (1-b-3*d)/4; (1-a-3*c)/4; (1-a-3*c)/4; 0; (1-a-3*c)/4; (1-b-3*d)/4; (1-b-3*d)/4; (1-a-3*c)/4; 0; (1-a-3*c)/4; (1-a-3*c)/4; (1-b-3*d)/4; (1-a-3*c)/4; 0; (1-b-3*d)/4; (1-b-3*d)/4; (1-b-3*d)/4; (1-b-3*d)/4];
        u1=[u1;u2];
    end
    u=u1
    f1=q./3+2*q1./3;
    f2=q./3+2*q2./3;
    f3=q./3+2*q3./3;
    f4=q./3;
    f5=q./3+q2./3;
    q=[f1;f2;f3;f4;f5];
end
toc;
figure;
xa=real(q);
ya=imag(q);
for i=1:4:length(q)
plot3([xa(i),xa(i+1),xa(i+2),xa(i+3),xa(i)],[ya(i),ya(i+1),ya(i+2),ya(i+3),ya(i)],[u(i),u(i+1),u(i+2),u(i+3),u(i)]);
hold on;
end
