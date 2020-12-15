
function gererate_result()
    data = circle(150, 0.1);
    [X,Y] = meshgrid(-6:1:6);
    data_move = [X(:)';Y(:)'];
    %data_move = circle(30, 1);
    
    %data_move = [(rand(1,10)-0.5)*10;(rand(1,10)-0.5)*12];
    
    subplot(1,3,1)
    %plot(data_move(1,:), data_move(2,:),'d','MarkerFaceColor', 'b')
    hold on
    ridge2(data, data_move, 0);
    axis([-6.5 6.5 -6.5 6.5])
    xlabel('KDE')
    
    subplot(1,3,2)
    %plot(data_move(1,:), data_move(2,:),'d','MarkerFaceColor', 'b')
    hold on
    ridge2(data, data_move, 1);
    axis([-6.5 6.5 -6.5 6.5])
    xlabel('KDE with LOG')
    
    subplot(1,3,3)
    %plot(data_move(1,:), data_move(2,:),'d','MarkerFaceColor', 'b')
    hold on
    ridge2(data, data_move, 2);
    axis([-6.5 6.5 -6.5 6.5])
    xlabel('LKDE with LOG')
end


function ridge2(data, data_move, algo)
    
    sigma = 0.9;
    epsion = 1e-6;
    max_iter = 3;
    
    for i = 1:size(data_move,2)
        track = [];
        for k = 1:max_iter
            if algo == 0
                [direction,temp1,temp2] = projection_gradient(data_move(:,i), data, sigma);
            elseif algo == 1
                [direction,temp1,temp2] = Hc(data_move(:,i),sigma,data,1);
            elseif algo == 2
                [direction,temp1,temp2] = Hc2(data_move(:,i),sigma,data,1);
            end
            track = [track, [data_move(:,i); angle(temp1, temp2);temp1;temp2]];
            data_move(:,i) = data_move(:,i)+direction;
            if norm(direction) < epsion
                break;
            end
        end
        %plot(track(1,:), track(2,:),'r-','Linewidth',1)
        %text(track(1,1)+0.2,track(2,1),num2str(track(3,1),'%1.2f'))
        quiver(track(1,1),track(2,1),0.5*track(4,1),0.5*track(5,1),'r')
        
        %annotation('arrow',[track(1,1) track(2,1)],[track(4,1) track(5,1)])
        hold on
        fprintf('iteration i=%d, iter=%d, error=%f\n', i, k, norm(direction));
    end
    
    plot(data(1,:),data(2,:),'bo')
    %hold on

    %plot(data_move(1,:),data_move(2,:),'b*')
    
end


function d = angle(a,b) 
    c = abs(a'*b/norm(a)/norm(b));
    d = sqrt(1-c^2);
end


function [g,temp1, temp2] = Hc(x, sigma, data, d)
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
    temp1 = V(:,2);
    temp2 = (c-x);
    %[U,~,~] = svd(H);
    %test_U = V(:,1:d)'*U;
    %[~,location] = min(sum(test_U.^2, 1));
    %P = U(:,location)*U(:,location)';
    P = eye(size(x, 1))-V(:,1:d)*V(:,1:d)';
    %g = 2*sum_r/size(data,2)/sigma^2*P*(c - x);
    g = 0.5*P*(c-x);
end


function [g, temp1, temp2] = Hc2(x, sigma, data, d)

    sq_distance = sum((data - x).^2,1);
    [~, ind] = sort(sq_distance);
    s = 3;
    
    H = zeros(size(x,1));
    c = zeros(size(x));
    sum_r = 0;
   % for i = 1:size(data,2)
    for i = ind(1:s)
        r = exp(-norm(x - data(:,i)).^2/(sigma^2));
       % H = H + r*(data(:,i)-x)*(data(:,i)-x)';
        c = c + r*data(:,i);
        sum_r = sum_r + r;
    end
    c = c/sum_r;
    B = zeros(size(x,1));
  %  for i = 1:size(data,2)
    for i = ind(1:s)
        r = exp(-norm(x - data(:,i)).^2/(sigma^2));
        B = B + r*(data(:,i)-c)*(data(:,i)-c)';
    end
    [V,~,~] = svd(B);
    temp1 = V(:,2);
    temp2 = (c-x);
    %[U,~,~] = svd(H);
    %test_U = V(:,1:d)'*U;
    %[~,location] = min(sum(test_U.^2, 1));
    %P = U(:,location)*U(:,location)';
    P = eye(size(x, 1))-V(:,1:d)*V(:,1:d)';
    %g = 2*sum_r/size(data,2)/sigma^2*P*(c - x);
    g = 0.5*P*(c-x);
end


function [P,temp1, g] = projection_gradient(x, Xi, sigma)
    g = gradient(x, Xi, sigma);
    H = Hession(x, Xi, sigma);
    %[U,D] = eig(H);
    [U,~,~] = svd(H);
    %[~, ind] = sort(diag(D),'descend');
    temp1 = U(:,2);
   % P = U(:,ind(2))*U(:,ind(2))'*g;
    P = U(:,2)*U(:,2)'*g;
end


function g = gradient(x, Xi, sigma)
    g = zeros(size(x));
    w = 0;
    for i = 1:size(Xi,2)
 %       g = g - 2./(sigma.^2)*exp(-norm(x-Xi(:,i)).^2/(sigma^2))*(x-Xi(:,i));
        g = g + exp(-norm(x-Xi(:,i)).^2/(sigma^2))*Xi(:,i);
        w = w + exp(-norm(x-Xi(:,i)).^2/(sigma^2));
    end
    g = g/w - x;
%    g = g / size(Xi,2); 
end


function H = Hession(x, Xi, sigma)
    s = size(x,1);
    H = zeros(s);
    for i = 1:size(Xi, 2)
  %      H = H + 2./(sigma.^4)*exp(-norm(x-Xi(:,i)).^2/(sigma^2))*...
  %          (2*(x-Xi(:,i))*(x-Xi(:,i))'-sigma.^2*eye(s));
        H = H + exp(-norm(x-Xi(:,i)).^2/(sigma^2))* (x-Xi(:,i))*(x-Xi(:,i))';
    end
 %   H = H / size(Xi,2);
end