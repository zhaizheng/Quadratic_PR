D = 3;
NumSample = 1000;
NumIni = 1000;
R = 2/3; r = 1/3;
sigma = 0.04;
RepeatTimes=5;
result = zeros(1,RepeatTimes);

for i = 1:RepeatTimes
    samples = torusUnif(NumSample, R, r);
    X2 = samples + sigma*randn(D,NumSample);
    data_ini = torusUnif(NumIni, R, r);
    X1 = data_ini+0.5*sqrt(sigma)/sqrt(D)*(2*rand(D,NumIni)-1);
    X = data_ini;
    [~, data4] = algorithm(X2, X1, 0.14);
    [P, d] = pdtorus(R, r, data4);
    result(i) = max_distance(data4, P);
    fprintf('Repeat %d, average error=%f:\n', i, mean(result(1:i)))
end


function max_d = max_distance(data1, data2)
    max_d = max(sqrt(sum((data1-data2).^2, 1)));
end


function [steps, data_move] = algorithm(data, data_move, sigma)
    epsion = 1e-6;
    max_iter = 10000;
    steps = 0;
    n = size(data_move,2);
    for i = 1:n
        for k = 1:max_iter
            direction = Hc5(data_move(:,i), sigma, data, 2);
            data_move(:,i) = data_move(:,i)+ direction;
            if norm(direction) < epsion
                break;
            end
        end
        steps = steps+k/n;
        %plot(track(1,:), track(2,:),'r-','Linewidth',1)
        %hold on
        %fprintf('iteration i=%d, iter=%d, error=%f\n', i, k, norm(direction));
    end
end


function g = Hc5(x,  sigma, data, d)
    c1 = C(data, sigma, x);
    c2 = C(data, sigma, c1);
    c3 = C(data, sigma, c2);
%     c4 = C(data, sigma, c3);
%     n = size(x,1);
%     E = eye(n) - (c3-c2)*(c4-c3)'/norm(c4-c3)^2;
%     [Q,~,~] = svd(E);
    B1 = zeros(size(x,1));
    B2 = zeros(size(x,1));
    for i = 1:size(data,2)
        %weight = (0.1+norm(x-c))/(0.2+norm(x-c));
        r1 = exp(-norm(c1 - data(:,i)).^2/(sigma^2));
        B1 = B1 + r1*(data(:,i)-c1)*(data(:,i)-c1)';
        
        r2 = exp(-norm(c2 - data(:,i)).^2/(sigma^2));
        B2 = B2 + r2*(data(:,i)-c2)*(data(:,i)-c2)';
    end
%     [W, D] = eig(Q(:,1:n-1)'*B*Q(:,1:n-1));
%     [~,ind] = sort(diag(D));
%     P = Q(:,1:n-1)*W(:,ind(1:n-d));
%     
    [V,~,~] = svd(B1+B2);
    P = eye(size(x, 1))-V(:,1:d)*V(:,1:d)';
    %P = eye(size(x, 1))-V(:,1:d)*V(:,1:d)';
    %g = 2*sum_r/size(data,2)/sigma^2*P*(c - x);
    g = 0.05*P*(c1 - x);
end


function c = C(data, sigma, x)
    c = zeros(size(x));
    sum_r = 0;
    for i = 1:size(data,2)
        r = exp(-norm(x - data(:,i)).^2/(sigma^2));
        c = c + r*data(:,i);
        sum_r = sum_r + r;
    end
    c = c/sum_r;
end