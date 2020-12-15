
function generate_result_digit()
    [data, data_move] = process_data();
    new_data = ridge2(data, data_move, 0);
    display(new_data, 7);
%     data = data_generation(0.02, 0.08);
%     data_move = data_generation(0.1, 0.3);
%     subplot(1,2,1)
%     plot(data_move(1,:),data_move(2,:),'d');
%     hold on
%     ridge2(data, data_move, 0)
%     subplot(1,2,2)
%     plot(data_move(1,:),data_move(2,:),'d');
%     hold on
%     ridge2(data, data_move, 1)
end

function [data, data_move] = process_data()
    [images,labels] = read_file();
    interest = images(:,:,labels==8);
    data = zeros(28*28,500);
    data_move = zeros(28*28,49);
    for i = 1:500
        data(:,i) = reshape(interest(:,:,i),[28*28,1]);
    end
    for i = 1:49
        data_move(:,i) = reshape(interest(:,:,i+500),[28*28,1]);
    end
end

function display(interest, num)
    I = zeros(28*num,28*num);
    for i = 1:num
        for j = 1:num
            I((i-1)*28+1:i*28, (j-1)*28+1:j*28) = interest(:,:,(i-1)*10+j);
        end
    end
    imshow(I)
end

function data_move = ridge2(data, data_move, algo)
    
    sigma = 0.25;
    epsion = 1e-5;
    max_iter = 10000;
    
    for i = 1:size(data_move,2)
        track = [];
        for k = 1:max_iter
            if algo == 1
                direction = projection_gradient(data_move(:,i), data, sigma);
            else
                direction = Hc(data_move(:,i),sigma,data,1);
            end
            data_move(:,i) = data_move(:,i)+0.01*direction;
            track = [track, data_move(:,i)];
            if norm(direction) < epsion
                break;
            end
        end
%         plot(track(1,:), track(2,:),'r-','Linewidth',1)
%         hold on
        fprintf('iteration i=%d, iter=%d, error=%f\n', i, k, norm(direction));
     end
%     
%     plot(data(1,:),data(2,:),'o')
%     hold on
% 
%     plot(data_move(1,:),data_move(2,:),'b*')
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