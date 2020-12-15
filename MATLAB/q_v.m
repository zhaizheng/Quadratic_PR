function q_value=q_v(X,Y,d,x)


N=length(X)-1;
 

for i=1:N
    SY(i)=Y(i+1)-Y(i);
end


Sq=0;
for i=1:N
    Sq=Sq+(Y(i)+Y(i+1)-2*d*Y(1))/(2*N);
end

 Move_Y=Y-Sq/(1-d);   %平移,使q积分值为0；
 

 k=floor((x-floor(x))*N)+1;
 q_value=(N.*(x-floor(x))-(k-1)).*SY(k)+Move_Y(k)-d*Move_Y(1);
end
