function f = compare1()
    [X, X1] = sphere();
    data1 = ridge2(X1, X1);
    fprintf('average distance error=%f\n',average_distance(data1, X));
    

end

function ave = average_distance(data1, data2)
    n = size(data1,2);
    sum_d = sum(sqrt(sum((data1-data2).^2, 1)),2);
    ave = sum_d/n;
end

function data_move = ridge2(data, data_move)
    

    epsion = 1e-5;
    max_iter = 10000;
    
    for i = 1:size(data_move,2)
        %track = [];
        for k = 1:max_iter
            direction = xia(data_move(:,i), data, 0.4, 3, 2);
            data_move(:,i) = data_move(:,i)+ direction;
            if norm(direction) < epsion
                break;
            end
        end
        %plot(track(1,:), track(2,:),'r-','Linewidth',1)
        %hold on
        fprintf('iteration i=%d, iter=%d, error=%f\n', i, k, norm(direction));
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
    direction = 0.1*P*(cx - x );
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