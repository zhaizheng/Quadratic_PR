
function run()
    repeat_num = 2;
    Sigma = [0.1, 0.2, 0.3];
    for j = 1:length(Sigma)
        sigma = Sigma(j);
        result = zeros([repeat_num, 6]);
        Steps = zeros([repeat_num, 3]);
        result2 = zeros([repeat_num, 6]);
        for i = 1:repeat_num
            [X, X1] = sphere2(sigma*sigma, 300);
            %show(X1, X);
            [steps1, data1] = ridge2(X1, X1, sigma, 0);
            [steps2, data2] = ridge2(X1, X1, 1.2*sigma, 1);
            [steps3, data3] = ridge2(X1, X1, 2*sigma, 2);
            result(i,:) = [average_distance(data1, X),...
                           average_distance(data2, X),...
                           average_distance(data3, X),...
                           average_distance(data1, bsxfun(@rdivide, data1, sqrt(sum(data1.^2,1)))),...
                           average_distance(data2, bsxfun(@rdivide, data2, sqrt(sum(data2.^2,1)))),...
                           average_distance(data3, bsxfun(@rdivide, data3, sqrt(sum(data3.^2,1)))),...
                           ];
            result2(i,:) = [max_distance(data1, X),...
                           max_distance(data2, X),...
                           max_distance(data3, X),...
                           max_distance(data1, bsxfun(@rdivide, data1, sqrt(sum(data1.^2,1)))),...
                           max_distance(data2, bsxfun(@rdivide, data2, sqrt(sum(data2.^2,1)))),...
                           max_distance(data3, bsxfun(@rdivide, data3, sqrt(sum(data3.^2,1)))),...
                           ];
            Steps(i,:) = [steps1,steps2,steps3];
        end
        fprintf('new average distance sigma=%f, error=%f, error2=%f, steps=%f\n',sigma, mean(result(:,1)),mean(result(:,4)), mean(Steps(:,1)));
        fprintf('KDE average distance sigma=%f, error=%f, error2=%f, steps=%f\n',sigma, mean(result(:,2)),mean(result(:,5)), mean(Steps(:,2)));
        fprintf('xia average distance sigma=%f, error=%f, error2=%f, steps=%f\n',sigma, mean(result(:,3)),mean(result(:,6)), mean(Steps(:,3)));
        fprintf('new max distance sigma=%f, error=%f, error2=%f, steps=%f\n',sigma, mean(result2(:,1)),mean(result2(:,4)), mean(Steps(:,1)));
        fprintf('KDE max distance sigma=%f, error=%f, error2=%f, steps=%f\n',sigma, mean(result2(:,2)),mean(result2(:,5)), mean(Steps(:,2)));
        fprintf('xia max distance sigma=%f, error=%f, error2=%f, steps=%f\n',sigma, mean(result2(:,3)),mean(result2(:,6)), mean(Steps(:,3)));
    end

end

function show(data1, X1)
    subplot(1,3,1)
    plot(data1(1,:),data1(2,:),data1(3,:),'o');
    hold on
    plot(X1(1,:), X1(2,:), X1(3,:),'*');
end


function ave = average_distance(data1, data2)
    n = size(data1,2);
    sum_d = sum(sqrt(sum((data1-data2).^2, 1)),2);
    ave = sum_d/n;
end

function max_d = max_distance(data1, data2)
    max_d = max(sqrt(sum((data1-data2).^2, 1)));
end

function [steps,data_move] = ridge2(data, data_move, sigma, algo)
    epsion = 1e-5;
    max_iter = 10000;
    steps = 0;
    n = size(data_move,2);
    for i = 1:n
        %track = [];
        for k = 1:max_iter
            if algo == 1
                direction = projection_gradient(data_move(:,i), data, sigma, 2);
            elseif algo == 0
                direction = Hc(data_move(:,i), sigma, data, 2);
            else  
                direction = xia(data_move(:,i), data, sigma, 3, 2);
            end
            data_move(:,i) = data_move(:,i)+ direction;
            %track = [track, data_move(:,i)];
            if norm(direction) < epsion
                break;
            end
        end
        steps = steps+k/n;
        %plot(track(1,:), track(2,:),'r-','Linewidth',1)
        %hold on
        %fprintf('iteration i=%d, iter=%d, error=%f', i, k, norm(direction));
    end
    
%     plot(data(1,:),data(2,:),'o')
%     hold on
% 
%     plot(data_move(1,:),data_move(2,:),'b*')
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
    n = size(data, 2);
    d = size(data, 1);
    indicator = ones([1,n]);
    indicator(ds>r^2) = 0;
    select = (data-data(:,i)).*indicator;
    cor = select*select';
    [U,~,~] = svd(cor);
    P = eye(d)-U(:,1:dim)*U(:,1:dim)';     
end

function g = Hc(x, sigma, data, d)
    c = zeros(size(x));
    sum_r = 0;
    for i = 1:size(data,2)
        r = exp(-norm(x - data(:,i)).^2/(sigma^2));
        c = c + r*data(:,i);
        sum_r = sum_r + r;
    end
    c = c/sum_r;
    B = zeros(size(x,1));
    %sigma = norm(x-c);
    for i = 1:size(data,2)
        r = exp(-norm(x - data(:,i)).^2/(sigma^2));
        B = B + r*(data(:,i)-c)*(data(:,i)-c)';
    end
    [V,~,~] = svd(B);
    P = eye(size(x, 1))-V(:,1:d)*V(:,1:d)';
    %g = 2*sum_r/size(data,2)/sigma^2*P*(c - x);
    g = 0.05*P*(c - x);
end

function P = projection_gradient(x, Xi, sigma, d)
    g = gradient(x, Xi, sigma);
    H = Hession(x, Xi, sigma);
    [U,~,~] = svd(H);
    P = U(:,1:d)*U(:,1:d)'*g;
end


function g = gradient(x, Xi, sigma)
    g = zeros(size(x));
    for i = 1:size(Xi,2)
        g = g - 2./(sigma.^2)*exp(-norm(x-Xi(:,i)).^2/(sigma^2))*(x-Xi(:,i));
    end
    g = g / size(Xi,2); 
end


function H = Hession(x, Xi, sigma)
    s = size(x,1);
    H = zeros(s);
    for i = 1:size(Xi, 2)
        H = H + 2./(sigma.^4)*exp(-norm(x-Xi(:,i)).^2/(sigma^2))*...
            (2*(x-Xi(:,i))*(x-Xi(:,i))'-sigma.^2*eye(s));
    end
    H = H / size(Xi,2);
end