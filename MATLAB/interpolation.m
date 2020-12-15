clear;
d = 2/3;
X = linspace(0, 1, 6);
Y=[2,100,-5,78,9,2];
%Y = randn(1,10);
N=length(X)-1;
outer_iteration = 10000;
inner_iteration = 300;
s = 0;
for i = 1:outer_iteration
    u = (i-1)/outer_iteration;
    si = 0;
    for j = 1:inner_iteration
        [~, v1, ~] = interpolation2(X, Y,d, N^j*u);
        si = si + d^j*v1;
    end
    [~, v2 ,~] = interpolation2(X, Y,d, u);
    s = s + si * v2/outer_iteration;
end
fprintf('int value:%E\n',s)





    x = 1.3;
    [xn, z, Move_Y] =  interpolation2(X, Y,d, x);
    plot(X, Move_Y-d*Move_Y(1),'-')
    hold on
    plot(xn, z,'o')


function [xn ,z, Move_Y] = interpolation2(X, Y,d, x)
    N = length(X)-1;
    square = 0;
    for i = 1:N
        square = square+(Y(i)+Y(i+1)-2*d*Y(1))/(2*N);
    end
    Move_Y = Y - square;
    find_k = floor(N*(x-floor(x)))+1;
%     lambda = (x-floor(x) - X(find_k+1))*(N);
    xn = x-floor(x);
%     z = lambda*(Move_Y(find_k+2)-Move_Y(find_k+1))+Move_Y(find_k+1);
z=(N*x-(find_k-1)).*(Move_Y(find_k+1)-Move_Y(find_k))+Move_Y(find_k)-d*Move_Y(1);
end