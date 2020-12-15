function ridge()
    data = data_generation(0.05, 0.3);
    sigma = 0.5;
    epsion = 1e-4;
    max_iter = 300;
    data_move = data_generation(0.1, 1);
    track = []
    for i = 1:size(data_move,2)
        for k = 1:max_iter
            direction = projection_gradient(data_move(:,i), data, sigma);
            data_move(:,i) = data_move(:,i)+direction;
            track = [track, data_move(:,i)];
            if norm(direction) < epsion
                break;
            end
        end
        
        fprintf('iteration i=%d, iter=%d, error=%f\n', i, k, norm(direction));
    end
    
    plot(data(1,:),data(2,:),'o')
    hold on
    plot(track(1,:), track(2,:),'.')
    %hold on
    %plot(data_move(1,:),data_move(2,:),'.')

end

function P = projection_gradient(x, Xi, sigma)
    g = gradient(x, Xi, sigma);
    H = Hession(x, Xi, sigma);
    [U,~,~] = svd(H);
    P = U(:,1)*U(:,1)'*g;
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