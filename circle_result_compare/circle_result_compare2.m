figure(1);
sample_size = 300; % in the paper, we let it to be 300;
D = 2; d = 1;
t = tiledlayout(2,2,'TileSpacing','Compact');
step1(sample_size, D, d);
result = step2(sample_size, D, d);

function step1(sample_size, D, d)
    [X1, X2, ~] = generate_sphere(0.08, D, sample_size, sample_size);
    neig = [100, 100, 100, 30, 30, 30];
    q = [0 -5 -10 0 -5 -10];
    sigma = 0.4;
    data_case = cell(1,6);
    for i = 1:6
        [~, data_case{i}]  =  algorithm(X2, X1, sigma, 0, d, neig(i), q(i));
    end
      
    for k = 1:2
        i = 2+(k-1)*3;
        nexttile
        plot(X2(1,:), X2(2,:), 'd');
        N2 = normalize(data_case{i});
        hold on
        plot(data_case{i}(1,:), data_case{i}(2,:), 'b.');
        hold on
        plot(N2(1,:), N2(2,:), 'r.');
        title(['{\it{l}}-SCRE: neig=',num2str(neig(i)),',q=',num2str(q(i))]);%d',neig(i),q(i));
        set(gca,'FontSize',16);
    end
end


function result = step2(sample_size, D, d)
    [X1, X2, ~] = generate_sphere(0.04, D, sample_size, sample_size);
    sigma = 0.1:0.1:0.9;
    result = zeros(8, length(sigma));
    data2 = cell(1,3);
    for i = 1:length(sigma)
        neig = 30; 
        q = 0;
        [~, data2{1}]  =  algorithm(X2, X1, sigma(i), 0, d, neig, q);
        q = -5;
        [~, data2{2}]  =  algorithm(X2, X1, sigma(i), 0, d, neig, q);
         q = -10;
        [~, data2{3}]  =  algorithm(X2, X1, sigma(i), 0, d, neig, q);
        [~, data6nn]   = algorithm(X2, X1, sigma(i), 8, d);
        [~, data6nnn]  = algorithm(X2, X1, sigma(i), 9, d);
        result(1,i) = average_distance(data2{1}, bsxfun(@rdivide, data2{1}, sqrt(sum(data2{1}.^2,1))));
        result(2,i) = average_distance(data2{2}, bsxfun(@rdivide, data2{2}, sqrt(sum(data2{2}.^2,1))));
        result(3,i) = average_distance(data2{3}, bsxfun(@rdivide, data2{3}, sqrt(sum(data2{3}.^2,1))));
        result(4,i) = average_distance(data6nn, bsxfun(@rdivide, data6nn, sqrt(sum(data6nn.^2,1))));
        result(5,i) = average_distance(data6nnn, bsxfun(@rdivide, data6nnn, sqrt(sum(data6nnn.^2,1))));
        result(6,i) = max_distance(data2{1}, bsxfun(@rdivide, data2{1}, sqrt(sum(data2{1}.^2,1))));
        result(7,i) = max_distance(data2{2}, bsxfun(@rdivide, data2{2}, sqrt(sum(data2{2}.^2,1))));
        result(8,i) = max_distance(data2{3}, bsxfun(@rdivide, data2{3}, sqrt(sum(data2{3}.^2,1))));
        result(9,i) = max_distance(data6nn, bsxfun(@rdivide, data6nn, sqrt(sum(data6nn.^2,1))));
        result(10,i) = max_distance(data6nnn, bsxfun(@rdivide, data6nnn, sqrt(sum(data6nnn.^2,1))));
    end
%     plot(result(1,:),'-*')
%     hold on
%     plot(result(2,:),'-o')
    %subplot('position',[0.04 0.14 0.45 0.82])
    marker1 = {'-.s','-.*','-.o','-->','--d'};
    nexttile
    for k = 1:5
        semilogy(sigma,result(k,:),marker1{k},'linewidth',1.2); 
        hold on; 
    end
    %semilogy(sigma, result(2,:),'--o','linewidth',1.2); %hold on; semilogy(sigma, result(3,:),'--s');
    legend('{\it{l}}-SCRE:\Gamma(0)','{\it{l}}-SCRE:\Gamma(-5)','{\it{l}}-SCRE:\Gamma(-10)','MFIT-i', 'MFIT-ii')
    xlabel('h')
    title('Average Margin')
    set(gca,'FontSize',16);
    
    
    %subplot('position',[0.53 0.14 0.45 0.82])
    nexttile
    for k = 1:5
        semilogy(sigma,result(5+k,:),marker1{k},'linewidth',1.2); 
        hold on; 
    end  
    %semilogy(sigma, result(6,:),'--o','linewidth',1.2); %hold on; semilogy(sigma, result(6,:),'--s')
    legend('{\it{l}}-SCRE:\Gamma(0)','{\it{l}}-SCRE:\Gamma(-5)','{\it{l}}-SCRE:\Gamma(-10)','MFIT-i', 'MFIT-ii')
    xlabel('h')
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
    data_ini = X + 0.5*sqrt(sigma)/sqrt(D)*(2*rand(D, NumIni)-1);
end


function [steps, data_move] = algorithm(data, data_move, sigma, algo, d, neig, q)
    epsion = 1e-10;
    max_iter = 10000;
    steps = 0;
    n = size(data_move,2);
    step = 0.5;
    for i = 1:n
        for k = 1:max_iter
            if algo == 0
                direction = local_mean_shift_parametrized(data_move(:,i), sigma, data, d, step, neig, q);
            elseif algo == 1
                direction = Hc1(data_move(:,i), sigma, data, d, step);
            elseif algo == 2
                direction = Hc2(data_move(:,i), sigma, data, d, step);
            elseif algo == 3
                direction = Hc3(data_move(:,i), sigma, data, d, step);
            elseif algo == 4
                direction = Hc44(data_move(:,i), sigma, data, d, step);
            elseif algo == 5
                direction = Hc55(data_move(:,i), sigma, data, d, step);
            elseif algo == 6
                direction = Hc66(data_move(:,i), sigma, data, d, step);
            elseif algo == 7
                direction = Hc77(data_move(:,i), sigma, data, d, step);
                %direction = xia(data_move(:,i), data, sigma, 3, 2);     
            elseif algo == 8
                direction = xia(data_move(:,i), data, sigma, 2, d, step);
            elseif algo == 9
                direction = mfit(data_move(:,i), data, sigma, 2, d, step);
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
    if alpha_sum >0 
        alpha = alpha_tilde./alpha_sum;
    else
        alpha = alpha_tilde;
    end
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
    if alpha_sum >0 
        alpha = alpha_tilde./alpha_sum;
    else
        alpha = alpha_tilde;
    end
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
    direction = step*P*(cx - x);
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


function g = local_mean_shift_parametrized(x, sigma, data, d, step, neig, q)

    sq_distance = sum((data - x).^2,1);
    [~, ind] = sort(sq_distance);
    
    c = zeros(size(x));
    sum_r = 0;
    %for i = 1:size(data,2)
    s = neig;
    for i = ind(1:s)
        r = exp(-norm(x - data(:,i)).^2/(sigma^2));
        c = c + r*data(:,i);
        sum_r = sum_r + r;
    end
    c = c/sum_r;
    B = zeros(size(x,1));
    for i = ind(1:s)
        r = exp(-norm(x - data(:,i)).^2/(sigma^2));
        B = B + r*(data(:,i)-c)*(data(:,i)-c)';
    end
    BB = B/sum_r-((1-q)*(c-x)*(c-x)');
    [U,E] = eig(BB);
    [~,ind] =sort(diag(E),'descend');
    V = U(:,ind);
    P = eye(size(x, 1))-V(:,1:d)*V(:,1:d)';
    g = step*P*(c - x);
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