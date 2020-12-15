
d = 0.4;
N = 9;
X = linspace(0, 1, 10);
Y = randn(1,10);
outer_iteration = 3000;
inner_iteration = 200;
s = 0;
for i = 1:outer_iteration
    u = (i-1)/outer_iteration;
    si = 0;
    for j = 1:inner_iteration
        [~, v1, ~] = interpolation2(X, Y, N^j*u);
        si = si + d^j*v1;
    end
    [~, v2 ,~] = interpolation2(X, Y, u);
    s = s + si * v2/outer_iteration;
end
fprintf('int value:%E\n',s)




function test_interpolation()
    X = linspace(0, 1, 10);
    Y = randn(1,10);
    x = 100.7;
    [xn, z, Move_Y] =  interpolation2(X, Y, x);
    plot(X, Move_Y,'-')
    hold on
    plot(xn, z,'o')
end


function [xn ,z, Move_Y] = interpolation2(X, Y, x)
    N = length(X);
    square = 0;
    for i = 1:N-1
        square = square+(Y(i)+Y(i+1))/2/(N-1);
    end
    Move_Y = Y - square;
    find_k = floor((N-1)*(x-floor(x)));
    lambda = (x-floor(x) - X(find_k+1))*(N-1);
    xn = x-floor(x);
    z = lambda*(Move_Y(find_k+2)-Move_Y(find_k+1))+Move_Y(find_k+1);
end