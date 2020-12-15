figure(1);
sample_size = 200; % in the paper, we let it to be 300;
step1(sample_size);
%figure(2);
result = step2(sample_size);
%saveas(gcf,'result1','eps');
%saveas(gcf,'result2','eps');

function step1(sample_size)
    [X1, X2, ~] = generate_sphere(0.08, 2, sample_size, sample_size);
    sigma = 0.4;
    [~, data6] = algorithm(X2, X1, sigma, 1);
    [~, data6n] = algorithm(X2, X1, sigma, 7);
%    [~, data6nn] = algorithm(X2, X1, sigma, 8);
    %subplot('position',[0.04 0.1 0.45 0.82]);
    subplot(2,2,1)
    plot(X2(1,:), X2(2,:), 'd');
    N2 = normalize(data6);
    hold on
    plot(data6(1,:), data6(2,:), 'b.');
    hold on
    plot(N2(1,:), N2(2,:), 'r.');
    title('SCRE')
    set(gca,'FontSize',16);
    
    %subplot('position',[0.53 0.1 0.45 0.82])
    subplot(2,2,2)
    plot(X2(1,:), X2(2,:), 'd');
    N3 = normalize(data6n);
    hold on
    plot(data6n(1,:), data6n(2,:), 'b.');
    hold on
    plot(N3(1,:), N3(2,:), 'r.');
    title('{\it{l}}-SCRE')
    set(gca,'FontSize',16);
%     
%     subplot(1,3,3)
%     plot(X2(1,:), X2(2,:), 'd');
%     N3 = normalize(data6nn);
%     hold on
%     plot(data6nn(1,:), data6nn(2,:), 'b.');
%     hold on
%     plot(N3(1,:), N3(2,:), 'r.');
%     title('Xia')
    

%    average_distance(data6, bsxfun(@rdivide, data6, sqrt(sum(data6.^2,1))))
%    average_distance(data6n, bsxfun(@rdivide, data6n, sqrt(sum(data6n.^2,1))))
end


function result = step2(sample_size)
    [X1, X2, ~] = generate_sphere(0.04, 2, sample_size, sample_size);
    sigma = 0.1:0.1:0.9;
    result = zeros(2, 9);
    for i = 1:9
        [~, data6] = algorithm(X2, X1, sigma(i), 1);
        [~, data6n] = algorithm(X2, X1, sigma(i), 7);
        [~, data6nn] = algorithm(X2, X1, sigma(i), 8);
        [~, data6nnn] = algorithm(X2, X1, sigma(i), 9);
        result(1,i) = average_distance(data6, bsxfun(@rdivide, data6, sqrt(sum(data6.^2,1))));
        result(2,i) = average_distance(data6n, bsxfun(@rdivide, data6n, sqrt(sum(data6n.^2,1))));
        result(3,i) = average_distance(data6nn, bsxfun(@rdivide, data6nn, sqrt(sum(data6nn.^2,1))));
        result(4,i) = average_distance(data6nn, bsxfun(@rdivide, data6nnn, sqrt(sum(data6nnn.^2,1))));
        result(5,i) = max_distance(data6, bsxfun(@rdivide, data6, sqrt(sum(data6.^2,1))));
        result(6,i) = max_distance(data6n, bsxfun(@rdivide, data6n, sqrt(sum(data6n.^2,1))));
        result(7,i) = max_distance(data6nn, bsxfun(@rdivide, data6nn, sqrt(sum(data6nn.^2,1))));
        result(8,i) = max_distance(data6nn, bsxfun(@rdivide, data6nnn, sqrt(sum(data6nnn.^2,1))));
    end
%     plot(result(1,:),'-*')
%     hold on
%     plot(result(2,:),'-o')
    %subplot('position',[0.04 0.14 0.45 0.82])
    subplot(2,2,3)
    semilogy(sigma,result(1,:),'-*','linewidth',1.2); 
    hold on; 
    semilogy(sigma, result(2,:),'--o','linewidth',1.2); %hold on; semilogy(sigma, result(3,:),'--s');
    legend('SCRE','{\it{l}}-SCRE')
    xlabel('h')
%     ylabel('margin')
    title('Average Margin')
    set(gca,'FontSize',16);
    %subplot('position',[0.53 0.14 0.45 0.82])
    subplot(2,2,4)
    semilogy(sigma,result(4,:),'-*','linewidth',1.2); 
    hold on; 
    semilogy(sigma, result(5,:),'--o','linewidth',1.2); %hold on; semilogy(sigma, result(6,:),'--s')
    legend('SCRE','{\it{l}}-SCRE')
    xlabel('h')
%     ylabel('margin')
    title('Hausdorff')
    set(gca,'FontSize',16);
end

function out = normalize(in)
    out = in*diag(1./sqrt(sum(in.^2)));
end

function ave = average_distance(data1, data2)
    n = size(data1,2);
    sum_d = sum(sqrt(sum((data1-data2).^2, 1)),2);
    ave = sum_d/n;
end

function max_d = max_distance(data1, data2)
    max_d = max(sqrt(sum((data1-data2).^2, 1)));
end

function [data_ini, samples, X] = generate_sphere(sigma, D, NumSample, NumIni)
    samples = randn(D, NumSample);
    samples = samples*diag(1./sqrt(sum(samples.^2)))+ sigma*randn(D, NumSample);
    
    data_ini = randn(D, NumIni);
    X = data_ini*diag(1./sqrt(sum(data_ini.^2)));
    data_ini = X + 0.5*sqrt(sigma)/sqrt(D)*(2*rand(D,NumIni)-1);
end

function [steps, data_move] = algorithm(data, data_move, sigma, algo)
    epsion = 1e-10;
    max_iter = 10000;
    steps = 0;
    n = size(data_move,2);
    step = 0.1;
    for i = 1:n
        for k = 1:max_iter
            if algo == 1
                direction = Hc1(data_move(:,i), sigma, data, 1, step);
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
                direction = Hc77(data_move(:,i), sigma, data, 1, step);
                %direction = xia(data_move(:,i), data, sigma, 3, 2);     
            elseif algo == 8
                direction = xia(data_move(:,i), data, sigma, 3, 1, step);
            elseif algo == 9
                direction = mfit(data_move(:,i), data, sigma, 3, 1, step);
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

function direction = mfit(x, data, r, beta, dim, step)
    d = size(data, 1);
    n = size(data, 2);
    d2 = 1 - sum((data-x).^2, 1)./(r^2);
    d2(d2<0) = 0;
    alpha_tilde = d2.^beta;
    alpha_sum = sum(alpha_tilde(:));
    alpha = alpha_tilde./alpha_sum;
    %cx = sum(data.* alpha, 2);
    Ns = zeros(d);
    c_vec = zeros(d,1);
    for i = 1:n
        if d2(i)>0
            ns = normal_space(data, i, r, dim);
            Ns = Ns+ns*alpha(i);
            c_vec = c_vec+alpha(i)*ns*(data(:,i)-x);
        end    
    end
    [U,~,~] = svd(Ns);
    P = U(:,1:d-dim)*U(:,1:d-dim)';
    direction = step*P*c_vec;
end

function direction = xia(x, data, r, beta, dim, step)
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
    direction = step*P*(cx - x );
end

function P = normal_space(data, i, r, dim)
    ds = sum((data-data(:,i)).^2, 1);
    ds(ds > r^2) = 0;
    n = size(data, 2);
    d = size(data, 1);
    indicator = zeros([1,n]);
    indicator(ds>0) = 1;
    select = (data-data(:,i)).*indicator;
    cor = select*select';
    [U,~,~] = svd(cor);
    P = eye(d)-U(:,1:dim)*U(:,1:dim)';     
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


function g = Hc77(x, sigma, data, d, step)

    sq_distance = sum((data - x).^2,1);
    [~, ind] = sort(sq_distance);
    
    c = zeros(size(x));
    sum_r = 0;
    %for i = 1:size(data,2)
    s = 20;
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