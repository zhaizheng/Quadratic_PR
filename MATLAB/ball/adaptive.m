
function adaptive()
%     sigma = 0.2;
% [X1, X2, X] = generate_sphere(0.04, 3, 1000, 1000);
% [~, data3] = algorithm(X2, X1, sigma, 3);
%     X = randn(2, 300);
%     Xn = X * diag(1./sqrt(sum(X.^2, 1)));
%     data = Xn + 0.1*randn(2, 300);
%     data2 = Xn + 0.1*randn(2, 300);
% %     [~, outputs1] = algorithm(data, data2, 0.2, 6, 0);
% %     [~, outputs2] = algorithm(data, data2, 0.2, 6, );
% %     subplot(1,2,1)
% %     plot(data(1,:), data(2,:), 'b.')
% %     hold on
%     for i = 1:10
%         [~ , outputs{i}] = algorithm(data, data2, i*0.1, 6, 0, 30);
%         subplot(2,5,i)
%         scatter(outputs{i}(1,:), outputs{i}(2,:))
%     end
    run()
end


function run()
    sigma = 0.2;
    repeat_num = 1;

    for i = 1:repeat_num
%         [X, X1] = sphere2(sigma*sigma, 300);
%         X2 = X1;
        [X1, X2, X] = generate_sphere(0.04, 3, 1000, 1000);
        %[~, X2] = sphere2(sigma*sigma, 300);
        %scatter3(X1(1,:),X1(2,:),X1(3,:));
        [~, data1] = algorithm(X2, X1, sigma, 1);
        [~, data2] = algorithm(X2, X1, sigma, 2);
        [~, data3] = algorithm(X2, X1, sigma, 3);
        [~, data4] = algorithm(X2, X1, sigma, 6);
        [~, data5] = algorithm(X2, X1, sigma, 5);
        [~, data6] = algorithm(X2, X1, sigma, 7);
        
        result(i,:) = [average_distance(data1, X),...
                       average_distance(data2, X),...
                       average_distance(data3, X),...
                       average_distance(data4, X),...
                       average_distance(data5, X),...
                       average_distance(data6, X),...
                       average_distance(data1, bsxfun(@rdivide, data1, sqrt(sum(data1.^2,1)))),...
                       average_distance(data2, bsxfun(@rdivide, data2, sqrt(sum(data2.^2,1)))),...
                       average_distance(data3, bsxfun(@rdivide, data3, sqrt(sum(data3.^2,1)))),...
                       average_distance(data4, bsxfun(@rdivide, data4, sqrt(sum(data4.^2,1)))),...
                       average_distance(data5, bsxfun(@rdivide, data5, sqrt(sum(data5.^2,1)))),...
                       average_distance(data6, bsxfun(@rdivide, data6, sqrt(sum(data6.^2,1)))),...
                       average_distance(X1, X),...
                       average_distance(X1, bsxfun(@rdivide, X1, sqrt(sum(X1.^2,1)))),...
                       ];
        result2(i,:) = [max_distance(data1, X),...
                       max_distance(data2, X),...
                       max_distance(data3, X),...
                       max_distance(data4, X),...
                       max_distance(data5, X),...
                       max_distance(data6, X),...
                       max_distance(data1, bsxfun(@rdivide, data1, sqrt(sum(data1.^2,1)))),...
                       max_distance(data2, bsxfun(@rdivide, data2, sqrt(sum(data2.^2,1)))),...
                       max_distance(data3, bsxfun(@rdivide, data3, sqrt(sum(data3.^2,1)))),...
                       max_distance(data4, bsxfun(@rdivide, data4, sqrt(sum(data4.^2,1)))),...
                       max_distance(data5, bsxfun(@rdivide, data5, sqrt(sum(data5.^2,1)))),...
                       max_distance(data6, bsxfun(@rdivide, data6, sqrt(sum(data6.^2,1)))),...
                       max_distance(X1, X),...
                       max_distance(X1, bsxfun(@rdivide, X1, sqrt(sum(X1.^2,1)))),...
                       ];
    end
    
    fprintf('r, sigma:%f ori error:%f HC1 error:%f  HC2:%f, HC3:%f, HC4:%f, Hc5:%f, Hc6:%f\n',...
        sigma, mean(result(:,13)), mean(result(:,1)), mean(result(:,2)), mean(result(:,3)), mean(result(:,4)),mean(result(:,6)), mean(result(:,5)));
    fprintf('r, sigma:%f ori error:%f HC1 error:%f  HC2:%f, HC3:%f, HC4:%f, Hc5:%f, Hc6:%f\n',...
        sigma, mean(result2(:,13)), mean(result2(:,1)), mean(result2(:,2)), mean(result2(:,3)),mean(result2(:,4)), mean(result2(:,6)), mean(result2(:,5)));
    
    fprintf('p, sigma:%f ori error:%f HC1 error:%f  HC2:%f, HC3:%f, HC4:%f, Hc5:%f, Hc6:%f\n',...
        sigma, mean(result(:,14)), mean(result(:,7)), mean(result(:,8)), mean(result(:,9)), mean(result(:,10)),mean(result(:,12)), mean(result(:,11)));
    fprintf('p, sigma:%f ori error:%f HC1 error:%f  HC2:%f, HC3:%f, HC4:%f, Hc5:%f, Hc6:%f\n',...
        sigma, mean(result2(:,14)), mean(result2(:,7)), mean(result2(:,8)), mean(result2(:,9)),mean(result2(:,10)), mean(result2(:,12)), mean(result2(:,11)));
   
end


function [data_ini, samples, X] = generate_sphere(sigma, D, NumSample, NumIni)
    samples = randn(3, NumSample);
    samples = samples*diag(1./sqrt(sum(samples.^2)))+ sigma*randn(3, NumSample);
    
    data_ini = randn(3, NumIni);
    X = data_ini*diag(1./sqrt(sum(data_ini.^2)));
    data_ini = X + 0.5*sqrt(sigma)/sqrt(D)*(2*rand(3,NumIni)-1);
end


function ave = average_distance(data1, data2)
    n = size(data1,2);
    sum_d = sum(sqrt(sum((data1-data2).^2, 1)),2);
    ave = sum_d/n;
end

function max_d = max_distance(data1, data2)
    max_d = max(sqrt(sum((data1-data2).^2, 1)));
end


function [steps, data_move] = algorithm(data, data_move, sigma, algo)
    epsion = 1e-8;
    max_iter = 10000;
    steps = 0;
    n = size(data_move,2);
    step = 0.3;
    for i = 1:n
        for k = 1:max_iter
            if algo == 1
                direction = Hc1(data_move(:,i), sigma, data, 2, step);
            elseif algo == 2
                direction = Hc2(data_move(:,i), sigma, data, 2, step);
            elseif algo == 3
                direction = Hc3(data_move(:,i), sigma, data, 2, step);
            elseif algo == 4
                direction = Hc44(data_move(:,i), sigma, data, 2, step);
            elseif algo == 5
                direction = Hc55(data_move(:,i), sigma, data, 2, step);
            elseif algo == 6
                direction = Hc66(data_move(:,i), sigma, data, 2, step);
            elseif algo == 7
                direction = Hc77(data_move(:,i), sigma, data, 2, step);
                %direction = xia(data_move(:,i), data, sigma, 3, 2);     
            end
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


function g = Hc1(x, sigma, data, d, step)
    c = zeros(size(x));
    sum_r = 0;
    for i = 1:size(data,2)
        r = exp(-norm(x - data(:,i)).^2/(sigma^2));
        c = c + r*data(:,i);
        sum_r = sum_r + r;
    end
    c = c/sum_r;
    B = zeros(size(x,1));
    for i = 1:size(data,2)
        r = exp(-norm(x - data(:,i)).^2/(sigma^2));
        B = B + r*(data(:,i)-x)*(data(:,i)-x)';
    end
    [V,~,~] = svd(B);
    P = eye(size(x, 1))-V(:,1:d)*V(:,1:d)';
    %g = 2*sum_r/size(data,2)/sigma^2*P*(c - x);
    g = step*P*(c - x);
end

function g = Hc2(x, sigma, data, d, step)
    c = zeros(size(x));
    sum_r = 0;
    for i = 1:size(data,2)
        r = exp(-norm(x - data(:,i)).^2/(sigma^2));
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
    P = eye(size(x, 1))-V(:,1:d)*V(:,1:d)';
    %g = 2*sum_r/size(data,2)/sigma^2*P*(c - x);
    g = step*P*(c - x);
end

function g = Hc3(x, sigma, data, d, step)

%     sq_distance = sum((data - x).^2,1);
%     [~, ind] = sort(sq_distance);
    
    c = zeros(size(x));
    sum_r = 0;
    for i = 1:size(data,2)
     
    %for i = ind(1:s)
        r = exp(-norm(x - data(:,i)).^2/(sigma^2));
        c = c + r*data(:,i);
        sum_r = sum_r + r;
    end
    c = c/sum_r;
    B = zeros(size(x,1));
    for i = 1:size(data,2)
%    for i = ind(1:s)
        r = exp(-norm(c - data(:,i)).^2/(sigma^2));
        B = B + r*(data(:,i)-c)*(data(:,i)-c)';
    end
    [V,~,~] = svd(B);
    P = eye(size(x, 1))-V(:,1:d)*V(:,1:d)';
    %g = 2*sum_r/size(data,2)/sigma^2*P*(c - x);
    g = step*P*(c - x);
end

function g = Hc55(x, sigma, data, d, step)

    sq_distance = sum((data - x).^2,1);
    [~, ind] = sort(sq_distance);
    
    c = zeros(size(x));
    sum_r = 0;
    %for i = 1:size(data,2)
    s = 15;
    for i = ind(1:s)
        r = exp(-norm(x - data(:,i)).^2/(sigma^2));
        c = c + r*data(:,i);
        sum_r = sum_r + r;
    end
    c = c/sum_r;
    B = zeros(size(x,1));
%    for i = 1:size(data,2)
    for i = ind(1:s)
        r = exp(-norm(x - data(:,i)).^2/(sigma^2));
        B = B + r*(data(:,i)-(c+x)/2)*(data(:,i)-(c+x)/2)';
    end
    [V,~,~] = svd(B);
    P = eye(size(x, 1))-V(:,1:d)*V(:,1:d)';
    %g = 2*sum_r/size(data,2)/sigma^2*P*(c - x);
    g = step*P*(c - x);
end

function g = Hc66(x, sigma, data, d, step)

    sq_distance = sum((data - x).^2,1);
    [~, ind] = sort(sq_distance);
    
    c = zeros(size(x));
    sum_r = 0;
    %for i = 1:size(data,2)
    s = 15;
    for i = ind(1:s)
        r = exp(-norm(x - data(:,i)).^2/(sigma^2));
        c = c + r*data(:,i);
        sum_r = sum_r + r;
    end
    c = c/sum_r;
    B = zeros(size(x,1));
%    for i = 1:size(data,2)
    for i = ind(1:s)
        r = exp(-norm(x - data(:,i)).^2/(sigma^2));
        B = B + r*(data(:,i)-x)*(data(:,i)-x)';
    end
    [V,~,~] = svd(B);
    P = eye(size(x, 1))-V(:,1:d)*V(:,1:d)';
    %g = 2*sum_r/size(data,2)/sigma^2*P*(c - x);
    g = step*P*(c - x);
end

function g = Hc77(x, sigma, data, d, step)

    sq_distance = sum((data - x).^2,1);
    [~, ind] = sort(sq_distance);
    
    c = zeros(size(x));
    sum_r = 0;
    %for i = 1:size(data,2)
    s = 15;
    for i = ind(1:s)
        r = exp(-norm(x - data(:,i)).^2/(sigma^2));
        c = c + r*data(:,i);
        sum_r = sum_r + r;
    end
    c = c/sum_r;
    B = zeros(size(x,1));
%    for i = 1:size(data,2)
    for i = ind(1:s)
        r = exp(-norm(x - data(:,i)).^2/(sigma^2));
        B = B + r*(data(:,i)-c)*(data(:,i)-c)';
    end
    [V,~,~] = svd(B);
    P = eye(size(x, 1))-V(:,1:d)*V(:,1:d)';
    %g = 2*sum_r/size(data,2)/sigma^2*P*(c - x);
    g = step*P*(c - x);
end


function g = Hc44(x, sigma, data, d, step)

%     sq_distance = sum((data - x).^2,1);
%     [~, ind] = sort(sq_distance);
    
    c = zeros(size(x));
    sum_r = 0;
    for i = 1:size(data,2)
        r = exp(-norm(x - data(:,i)).^2/(sigma^2));
        c = c + r*data(:,i);
        sum_r = sum_r + r;
    end
    c = c/sum_r;
    B1 = zeros(size(x,1));
    B2 = zeros(size(x,1));
    for i = 1:size(data,2)
    %for i = ind(1:s)
        r1 = exp(-norm(c - data(:,i)).^2/(sigma^2));
        r2 = exp(-norm(x - data(:,i)).^2/(sigma^2));
        B1 = B1 + r1*(data(:,i)-c)*(data(:,i)-c)';
        B2 = B2 + r2*(data(:,i)-x)*(data(:,i)-x)';
    end
    [V,~,~] = svd(B1+B2);
    P = eye(size(x, 1))-V(:,1:d)*V(:,1:d)';
    %g = 2*sum_r/size(data,2)/sigma^2*P*(c - x);
    g = step*P*(c - x);
end

function g = Hc33(x, sigma, data, d, step)
    c = zeros(size(x));
    n = size(data, 2);
    sum_r = 0;
    for i = 1:size(data,2)
        r = exp(-norm(x - data(:,i)).^2/(sigma^2));
        c = c + r*data(:,i);
        sum_r = sum_r + r;
    end
    c = c/sum_r;
    B1 = zeros(size(x,1));
    B2 = zeros(size(x,1));
    for i = 1:size(data,2)
        %weight = (0.1+norm(x-c))/(0.2+norm(x-c));
        %r1 = exp(-norm(c - data(:,i)).^2/(sigma^2));
        r2 = exp(-norm(x - data(:,i)).^2/(sigma^2));
        B1 = B1 + r2*(data(:,i)-c)*(data(:,i)-c)';
        B2 = B2 + r2*(data(:,i)-x)*(data(:,i)-x)';
    end
    [V1,~,~] = svd(B1); P1 = V1(:,1:d)*V1(:,1:d)';
    [V2,~,~] = svd(B2); P2 = V2(:,1:d)*V2(:,1:d)';
    mu11 = (c-x)'*P1*(c-x); mu12 = (c-x)'*(eye(size(x,1))-P1)*(c-x);
    mu21 = (c-x)'*P2*(c-x); mu22 = (c-x)'*(eye(size(x,1))-P2)*(c-x);
    if mu12<mu11 && mu22<mu21
        P = eye(size(x, 1))-P2;
    else
        P = eye(size(x, 1))-P1;
    end
    %g = 2*sum_r/size(data,2)/sigma^2*P*(c - x);
    g = step*P*(c - x);
end


function g = Hc6(x, sigma, data, d, step)
    c = zeros(size(x));
    n = size(data, 2);
    sum_r = 0;
%     squared = sum((data-x).*(data-x), 1);
%     [~,ind] = sort(squared);
    for i = 1:n
        r = exp(-norm(x - data(:,i)).^2/(sigma^2));
        c = c + r*data(:,i);
        sum_r = sum_r + r;
    end
    c = c/sum_r;
    B1 = zeros(size(x,1));
    %B2 = zeros(size(x,1));
    %m = (x+c)./2;
    for i = 1:n
        %weight = (0.1+norm(x-c))/(0.2+norm(x-c));
        %r1 = exp(-norm(c - data(:,i)).^2/(sigma^2));
        %r2 = exp(-norm(x - data(:,i)).^2/(sigma^2));
        r3 = exp(-norm(c - data(:,i)).^2/(sigma^2));
        %r4 = exp(-norm(m - data(:,i)).^2/(sigma^2));
        %B1 = B1 + r3*(data(:,i)-c)*(data(:,i)-c)'+r2*(data(:,i)-x)*(data(:,i)-x)'+r4*(data(:,i)-m)*(data(:,i)-m)';
        B1 = B1 + r3*(data(:,i)-c)*(data(:,i)-c)';
        %B2 = B2 + r3*(data(:,i)-c)*(data(:,i)-c)';
    end
    [V1,~,~] = svd(B1); P1 = V1(:,1:d)*V1(:,1:d)';
    P1 = eye(size(x, 1))-P1;
    cn = x+P1*(c-x);
    B2 = zeros(size(x,1));
    for i = 1:n
        r4 = exp(-norm(cn-data(:,i)).^2/(sigma^2));
        B2 = B2+r4*(data(:,i)-cn)*(data(:,i)-cn)';
    end
    
    [V2,~,~] = svd(B2); P2 = V2(:,1:d)*V2(:,1:d)';
    P2 = eye(size(x, 1))-P2;
    %g = 2*sum_r/size(data,2)/sigma^2*P*(c - x);
    g = step*P2*(cn-x);
    %g2 = 0.05*P2*(c - x);
    %g = (g1+g2)/2;
end

function g = Hc7(x, sigma, data, d, step)
    c = zeros(size(x));
    n = size(data, 2);
    sum_r = 0;
    for i = 1:n
        r = exp(-norm(x - data(:,i)).^2/(sigma^2));
        c = c + r*data(:,i);
        sum_r = sum_r + r;
    end
    c = c/sum_r;
    %%
    B1 = zeros(size(x,1));
    for i = 1:n
        r1 = exp(-norm(c - data(:,i)).^2/(sigma^2));
        B1 = B1 + r1*(data(:,i)-c)*(data(:,i)-c)';
    end
    [V1,~,~] = svd(B1); P1 = V1(:,1:d)*V1(:,1:d)';
    P1 = eye(size(x, 1))-P1;
    cn = x+P1*(c-x);
    %%
    B2 = zeros(size(x,1));
    for i = 1:n
        r2 = exp(-norm(cn-data(:,i)).^2/(sigma^2));
        B2 = B2+r2*(data(:,i)-cn)*(data(:,i)-cn)';
    end
    [V2,~,~] = svd(B2); P2 = V2(:,1:d)*V2(:,1:d)';
    P2 = eye(size(x, 1))-P2;
    cn2 = x+P2*(cn - x);
    %%
    B3 = zeros(size(x,1));
    for i = 1:n
        r3 = exp(-norm(cn2-data(:,i)).^2/(sigma^2));
        B3 = B3+r3*(data(:,i)-cn2)*(data(:,i)-cn2)';
    end
    [V3,~,~] = svd(B3); P3 = V3(:,1:d)*V3(:,1:d)';
    P3 = eye(size(x, 1))-P3;
    %%
    g = step*P3*(cn2-x);
end

function g = Hc4(x, sigma, data, d, step)
%     kk = 10;
    sq_distance = sum((data - x).^2,1);
    [~, ind] = sort(sq_distance);
%     sdata = (data(:,ind(1:kk*d))-x)'*(data(:,ind(1:kk*d))-x);
%     [U,L,~] = svd(sdata);
% %     IS = diag(L)>1e-6;
% %     para = U(:,IS)*diag(1./(diag(L(IS,IS))))*U(:,IS)'*ones(3*d,1);
%     para = U/(L+eye(size(L))*1e-3)*U'*ones(kk*d,1);
%     w = para./sum(para);
%     c = data(:, ind(1:kk*d))*w;
    c = data(:, ind(1));
    B = zeros(size(x,1));
    for i = 1:size(data,2)
        %weight = (0.1+norm(x-c))/(0.2+norm(x-c));
        r = exp(-norm(c - data(:,i)).^2/(sigma^2));
        B = B + r*(data(:,i)-c)*(data(:,i)-c)';
    end
    [V,~,~] = svd(B);
    P = eye(size(x, 1))-V(:,1:d)*V(:,1:d)';
    %g = 2*sum_r/size(data,2)/sigma^2*P*(c - x);
    g = step*P*(c - x);
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


function g = Hc5(x,  sigma, data, d, step)
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
    g = step*P*(c1 - x);
end

function direction = xia(x, data, r, beta, dim)
    d = size(data, 1);
    n = size(data, 2);
    d2 = 1 - sum((data-x).^2, 1)./(r^2);
    d2(d2<0) = 0;
    alpha_tilde = d2.^beta;
    alpha_sum = sum(alpha_tilde(:));
    alpha = alpha_tilde./alpha_sum;
    cx = sum(data.* alpha, 2);
    Ns = zeros(d);
    for i = 1:n
        if d2(i)>0
            ns = normal_space(data, i, r, dim);
            Ns = Ns+ns*alpha(i);
        end    
    end
    [U,~,~] = svd(Ns);
    P = U(:,1:d-dim)*U(:,1:d-dim)';
    direction = 0.05*P*(cx - x );
end

function P = normal_space(data, i, r, dim)
    ds = sum((data-data(:,i)).^2, 1);
    ds(ds > r^2) = 0;
    n = size(data, 2);
    d = size(data, 1);
    indicator = ones([1,n]);
    indicator(ds>r^2) = 0;
    select = (data-data(:,i)).*indicator;
    cor = select*select';
    [U,~,~] = svd(cor);
    P = eye(d)-U(:,1:dim)*U(:,1:dim)';     
end
