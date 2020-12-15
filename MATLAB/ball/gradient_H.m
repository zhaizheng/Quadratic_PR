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
            Rj = Rj + (-2*exp(-norm(x-data(:,i)).^2/sigma^2)*t/(sigma^2))*(x-data(:,i))*(x-data(:,i))'-sigma^2/2*eye(D,D);
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