clc,clear 
syms y 

% (1)% X=[0,1/3,2/3,1];
%    Y=[0,9,-5,0];
%    d=1/2;

X=[0,1/5,2/5,3/5,4/5,1];
  Y=[2,100,-5,78,9,2];
 d=2/3;


N=length(X)-1;


for i=1:N
    SY(i)=Y(i+1)-Y(i);
end


Sq=0;
for i=1:N
    q2(i)=(N.*y-(i-1)).*SY(i)+Y(i)-d*Y(1);
    Sq=Sq+int(q(i),(i-1)/N,i/N);
end

 Y=Y-eval(Sq)/(1-d);   %平移,使q积分值为0；
 

Sq=0;
for i=1:N
    q2(i)=(N.*y-(i-1)).*SY(i)+Y(i)-d*Y(1);
    Sq=Sq+int(q(i),(i-1)/N,i/N);
end
Sq  %输出q的积分

% (1) % q_1=@(t)((N.*t-(1-1)).*SY(1)+Y(1)-d*Y(1)).*(0<=t<1/N)+((N.*t-(2-1)).*SY(2)+Y(2)-d*Y(1)).*(1/N<=t<2/N)+((N.*t-(3-1)).*SY(3)+Y(3)-d*Y(1)).*(2/N<=t<=3/N);

q_1=@(t)((N.*t-(1-1)).*SY(1)+Y(1)-d*Y(1)).*(0<=t<1/N)+((N.*t-(2-1)).*SY(2)+Y(2)-d*Y(1)).*(1/N<=t<2/N)+((N.*t-(3-1)).*SY(3)+Y(3)-d*Y(1)).*(2/N<=t<3/N)+((N.*t-(4-1)).*SY(4)+Y(4)-d*Y(1)).*(3/N<=t<4/N)+((N.*t-(5-1)).*SY(5)+Y(5)-d*Y(1)).*(4/N<=t<=5/N);



%分割求和取极限求积分

I=0;
for m=1:10000
    sum=0;
    for i=1:300
        sum=sum+d^i*q(mod((N^i)*(m/10000),1));
    end
        I=I+(1/10000).*q(m/10000).*sum;
end
I   %输出左边积分值


% 
% I=0;
% for m=1:10000
%     sum=0;
%     for i=1:300
%         sum=sum+d^i*q_1(mod((N^i)*(m/10000),1));
%     end
%         I=I+(1/10000).*q_1(m/10000).*sum;
% end
% I   %输出左边积分值

II=0;
for i=1:N
    II=II+int(q(i).^2,(i-1)/N,i/N);
end

eval(II)  %输出右边积分值

if I<=eval(II)
    fprintf('Beautiful')
end



function q_value=q(x)
 k=floor(x*N)+1;
 q_value=(N.*x-(k-1)).*SY(k)+Y(k)-d*Y(1)
end