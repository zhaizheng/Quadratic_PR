t = -1:0.01:1;
t2 = -1.5:0.1:1.2;
y = sin(4*t)+0.08*randn(size(t));
plot(t, sin(4*t), '--', 'LineWidth', 1)
hold on
plot(t, y, 'o','MarkerSize', 4)
hold on
x = [-0.2; 0.4];
plot(x(1), x(2), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'b')
% t = 0:0.1:pi;
% x = 0.8*cos(t);
% y = 0.8*sin(t)+0.3;
% plot(x,y,'-','LineWidth',2)
%
%%
sigma = 0.3;
u = covariance_matrix(x, x, [t;y], sigma, 1);
T = x+u*t2;
hold on
plot(T(1,:),T(2,:), 'b-')
%%
m = shift_mean(x, [t;y], sigma);
u2 = covariance_matrix(m, m, [t;y], sigma, 1);
T2 = m+u2*t2;
hold on
plot(T2(1,:),T2(2,:), 'r*-')
hold on
plot(m(1),m(2), 'o','MarkerSize', 8, 'MarkerFaceColor', 'r')

%%
m = shift_mean(x, [t;y], sigma);
u3 = covariance_matrix2(m, m, [t;y], sigma, 1);
T3 = m+u3*t2;
hold on
plot(T3(1,:),T3(2,:), 'm.-')
hold on
plot(m(1),m(2), 'o','MarkerSize', 8, 'MarkerFaceColor', 'r')

axis([-1, 1,-1.5, 1.5])

%%
function u = covariance_matrix(x1, x2, data, sigma, s)
    d = size(data, 1);
    n = size(data, 2);
    C = zeros(d, d);
    for i = 1:n
        C = C + exp(-norm(data(:,i)-x1)^2/sigma^2)*(data(:,i)-x2)*(data(:,i)-x2)';
    end
    [U1,~,~] = svd(C);
    u = U1(:,1:s);
end


function u = covariance_matrix2(x1, x2, data, sigma, s)
    d = size(data, 1);
    n = size(data, 2);
    C = zeros(d, d);
    for i = 1:n
        C = C + exp(-norm(data(:,i)-x1)^2/sigma^2)*(data(:,i)-x2)*(data(:,i)-x2)'/(norm((data(:,i)-x2))^2);
    end
    [U1,~,~] = svd(C);
    u = U1(:,1:s);
end

function m = shift_mean(x, data, sigma)
    d = size(data, 1);
    n = size(data, 2);
    m = zeros(d, 1);
    w = 0;
    for i = 1:n
        m = m + exp(-norm(data(:,i)-x)^2/sigma^2)*data(:,i);
        w = w + exp(-norm(data(:,i)-x)^2/sigma^2);
    end
    m = m/w;
end
