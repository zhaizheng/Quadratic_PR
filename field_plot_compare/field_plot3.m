
demo = 1;% demo = 2
figure
t = tiledlayout(1,4,'TileSpacing','Compact');
generate_result(demo);



function generate_result(demo)
    data = circle(150, 0.1);
    [X,Y] = meshgrid(-6:1:6);
    data_move = [X(:)';Y(:)'];
    if demo == 2
        data_move = circle(30, 1);
    end
  
    Title = {'\Gamma(-1)','\Gamma(0)','\Gamma(0.5)','\Gamma(1)'};
    for i = 1:4
        %subplot('position', [0.04*i+0.29*(i-1) 0.12 0.28 0.86]);
        nexttile
        if demo==2
            plot(data_move(1,:),data_move(2,:),'d','MarkerFaceColor', 'r');           
        end
        ridge2(data, data_move, i-1, demo);
        axis([-6.5 6.5 -6.5 6.5])
        title(Title{i},'Interpreter','tex')
        set(gca,'FontSize',16);
    end
end


function ridge2(data, data_move, algo, demo)
    
    sigma = 0.9;
    epsion = 1e-6;
    max_iter = 3;
    step = 0.5;
    
    for i = 1:size(data_move,2)
        track = [];
        for k = 1:max_iter
            if algo == 0
                q = -1;
            elseif algo == 1
                q = 0;
            elseif algo == 2
                q = 0.5;
            elseif algo == 3
                q = 1;
            end
            [direction,temp1,temp2] = Hc(data_move(:,i),sigma,data,1,step, q);
            track = [track, [data_move(:,i); angle(temp1, temp2); temp1; temp2]];
            data_move(:,i) = data_move(:,i)+direction;
            if norm(direction) < epsion
                break;
            end
        end
        %plot(track(1,:), track(2,:),'r-','Linewidth',1)
        if demo==2
            text(track(1,1)+0.2,track(2,1),num2str(track(3,1),'%1.2f'), 'FontSize',14)
        else
            h1 = quiver(track(1,1),track(2,1),0.5*track(4,1),0.5*track(5,1),'linewidth',2);
            set(h1,'AutoScale','off', 'AutoScaleFactor', 2);
        end
        
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


function [g,temp1, temp2] = Hc(x, sigma, data, d, step, q)
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
        B = B + r*(data(:,i)-x)*(data(:,i)-x)';
    end
    BB = B/sum_r-((1-q)*(c-x)*(c-x)');
    [U,E] = eig(BB);
    [~,ind] =sort(diag(E),'descend');
    V = U(:,ind);
    
    temp1 = V(:,2);
    temp2 = (c-x);

    P = eye(size(x, 1))-V(:,1:d)*V(:,1:d)';
    g = step*P*(c-x);
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

function data = circle(n,sigma)
    %n = 1000 sigma = 0.2
    t = linspace(pi/2,5*pi,n);
    r = linspace(1,7,n);
    
    x = r.*cos(t)+sigma*randn(1,n);
    y = r.*sin(t)+sigma*randn(1,n);
    data = [x;y];
%     plot(x,y,'.')
end
