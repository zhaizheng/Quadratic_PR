function test_n()
    sigma = 0.2;
    [X, X1] = sphere2(sigma*sigma, 300);
    epsion = 1e-5;
    max_iter = 10000;
    steps = 0;
    n = size(X, 2);
    data_move = X1;
    d = 2;
    for i = 1:n
        for k = 1:max_iter
            [G1, G2] = gradient_H(data_move(:,i), X1, sigma, d);
            direction = -2*G1;
            data_move(:,i) = data_move(:,i)+ direction;
            if norm(direction) < epsion
                break;
            end
        end
        steps = steps+k/n;
        fprintf('samples:%d, steps:%d\n',i, k);
    end
end

function ave = average_distance(data1, data2)
    n = size(data1,2);
    sum_d = sum(sqrt(sum((data1-data2).^2, 1)),2);
    ave = sum_d/n;
end

function max_d = max_distance(data1, data2)
    max_d = max(sqrt(sum((data1-data2).^2, 1)));
end

function [G1, G2] = gradient_H(x, data, sigma, d)
    n = size(data, 2); D = size(x, 1);
    H = zeros(D, D);
    for i = 1:n
        H = H + exp(-norm(x-data(:,i)).^2/sigma^2)*((x-data(:,i))*(x-data(:,i))'-eye(D,D)*sigma^2/2);
    end
    [U, ei] = eig(H);
    [~, indx] = sort(diag(ei));
    U2 = U(:,indx(1:D-d));
    UU2 = U2*U2';
    PH_s = zeros(D, 1);
    Tr_s = zeros(D, 1);
    for j = 1:D
        Rj = zeros(D, D);
        for i = 1:n
            t = x(j)-data(j,i);
            Rj = Rj + (-2*exp(-norm(x-data(:,i)).^2/sigma^2)*t/(sigma^2))*((x-data(:,i))*(x-data(:,i))'-sigma^2/2*eye(D,D));
            T = zeros(D,D);
            T(j,:) = (x-data(:,i))';
            Rj = Rj + exp(-norm(x-data(:,i)).^2/sigma^2)* (T + T' + t*eye(D,D));
        end
        PH_s(j) = UU2(:)'*Rj(:);
        Tr_s(j) = trace(Rj);
    end
    G1 = PH_s;
    G2 = Tr_s;
end