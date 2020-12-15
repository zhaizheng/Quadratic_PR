function boundary()
    data = data_generation(0.05, 0.3);
    sigma = 4;
    epsion = 1e-6;
    max_iter = 300;
    data_move = data_generation(0.03, 1);
    data1 = data(:,1:size(data,2)/2);
    data2 = data(:,size(data,2)/2+1:size(data,2));
    for i = 1:size(data_move,2)
        for k = 1:max_iter
            [p1,~] = projection_gradient(data_move(:,i), data1, sigma);
            [p2,~] = projection_gradient(data_move(:,i), data2, sigma);
            u1 = (p1+p2)'*p1;
            u2 = (p1+p2)'*p2;
            if u1 > u2 
                data_move(:,i) = data_move(:,i)- (u1-u2)*p1;
            else
                data_move(:,i) = data_move(:,i)- (u2-u1)*p2;
            end
            
            %data_move(:,i) = data_move(:,i)+direction;
            if abs(u1-u2) < epsion
                break;
            end
        end
        fprintf('iteration i=%d, iter=%d, error=%f\n', i, k, abs(u1-u2));
    end
    plot(data(1,:),data(2,:),'o')
    hold on
    plot(data_move(1,:), data_move(2,:),'.')
    
end

function [point, P] = projection_gradient(x, Xi, sigma)
    g = gradient(x, Xi, sigma);
    H = Hession(x, Xi, sigma);
    [U,~,~] = svd(H);
    P = U(:,1)*U(:,1)';
    point = U(:,1)*U(:,1)'*g;
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