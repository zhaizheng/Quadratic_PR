
%[l,r] = find_interval(0.25, 100)
summation(linspace(1,3,5))

function [left, right] = find_interval(x, N)
    left = floor(x*N);
    right = left+1;
end

function s = summation(x)
    s = sum(x(:));
end

function test_li2()
    hahaha(1)
end

function q=hahaha(t)
    t=mod(t,1);
    q=q_1(t); 
end

function q = q_1(t)
        q=@(t)((N.*t-(1-1)).*SY(1)+Y(1)-d*Y(1)).*(0<=t<1/N)+((N.*t-(2-1)).*SY(2)+Y(2)-d*Y(1)).*(1/N<=t<2/N)+((N.*t-(3-1)).*SY(3)+Y(3)-d*Y(1)).*(2/N<=t<=3/N);
end


