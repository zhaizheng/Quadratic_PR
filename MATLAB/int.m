clear 
syms t 

X=[0,1/6,2/6,3/6,4/6,5/6,1];
   Y1 = randn(1,6);
Y=[Y1,Y1(1)]  %% 分形插值函数 f(0)=f(1)
 d=3/4;  %%%  纵向尺度因子 0<d<1, d*N>1


N=length(X)-1;

for i=1:N
    delta_Y(i)=Y(i+1)-Y(i);
end

%%%利用梯形面积计算q(t)在[0,1]上的积分值
Sq=0;
for i=1:N
    Sq=Sq+(Y(i)+Y(i+1)-2*d*Y(1))/(2*N);
end

Y=Y-Sq/(1-d);   %平移Y,使q积分值为0；
 

%%%平移Y后，计算q的积分值
Sq=0;
for i=1:N
    q(i)=(N.*t-(i-1)).*delta_Y(i)+Y(i)-d*Y(1); 
    Sq=Sq+int(q(i),(i-1)/N,i/N);
end
Sq  %验证q的积分是否为0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%q(Nu)与q(u)的乘积，在[0，1]上的积分是否为0.

III=0;
for i=1:50000
    III=III+(1/50000).*q_v(X,Y,d,N*i/50000).*q_v(X,Y,d,i/50000);
end
III

SS2=0;
for i=1:N
    SS1=0;
    for j=1:N
        Temper_q=t*((N^2*t-(i-1)*N-(j-1))*delta_Y(j)+Y(j)-d*Y(1));
        SS1=SS1+int(Temper_q,(i-1)/N+(j-1)/(N^2),(i-1)/N+j/(N^2));
    end
    SS2=SS2+N*delta_Y(i)*SS1;
end
T_SS2=eval(SS2);
fprintf('the interpolation of q(Nu)q(u) is%6.2f\n',T_SS2)

SS4=0;
for i=1:N
    SS3=0;
    for j=1:N
        SS3=SS3+int(q(j)*t,(j-1)/N,j/N);
    end
    SS4=SS4+delta_Y(i)*SS3/N;   
end
T_SS4=eval(SS4);
fprintf('the interpolation of q(Nu)q(u) is%6.2f\n',T_SS4)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%用黎曼和求左边积分

I=0;
for m=1:50000
    S1=0;
    for i=1:300
        S1=S1+d^i*q_v(X,Y,d,mod((N^i)*(m/50000),1));
    end
        I=I+(1/50000).*q_v(X,Y,d,m/50000).*S1;
end
I   %输出左边积分值



II=0;
for i=1:N
    II=II+int(q(i).^2,(i-1)/N,i/N);
end

eval(II)  %输出q^2的积分值


if I<=eval(II)
    fprintf('Beautiful')
end



%%%函数q(x)
function q_value=q_v(X,Y,d,x)

N=length(X)-1;
 
for i=1:N
    delta_Y(i)=Y(i+1)-Y(i);
end

 k=floor((x-floor(x))*N)+1;  %%% 找到x落在那个区间
 q_value=(N.*(x-floor(x))-(k-1)).*delta_Y(k)+Y(k)-d*Y(1); %%计算q(x)
end
