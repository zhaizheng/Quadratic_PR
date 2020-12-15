function ridge_find()

end


function g = Hc(x, sigma, data, d)
    H = zeros(size(x,1));
    c = zeros(size(x));
    sum_r = 0;
    for i = 1:size(data,2)
        r = exp(-norm(x - data(:,i)).^2/(sigma^2));
        H = H + r*(data(:,i)-x)*(data(:,i)-x)';
        c = c + r*data(:,i);
        sum_r = sum_r + r;
    end
    c = c/sum_r;
    B = zeros(size(x,1));
    for i = 1:size(data,2)
        r = exp(-norm(x - data(:,i)).^2/(sigma^2));
        B = B + r*(data(:,i)-c)*(data(:,i)-c)';
    end
    [V,~,~] = svd(B);
    [U,~,~] = svd(H);
    test_U = V(:,1:d)'*U;
    [~,location] = min(sum(test_U.^2, 1));
    P = U(:,location)*U(:,location)';
    g = 2*sum_r/size(data,2)/sigma^2*P*(c - x);
end