sample_size = 3000; % in the paper, we let it to be 300;
D = 2; d = 1;
%[X1, X2, T] = generate_random(D, sample_size, 300);
[X1, X2, T] = generate_sphere(0.4, D, sample_size, 300);


Re = cell(1,5);
sign = cell(1,5);
result = zeros(1,5);
t = tiledlayout(1,5,'TileSpacing','Compact');
Q = [-10000, -1, 0, 0.5, 1];
for j = 1:5
    [~, Re{j},sign{j}] = algorithm(X2, X1, 0.3, Q(j), d); %(j-1)*0.25
    nexttile
    plot(Re{j}(1,sign{j}<0.045),Re{j}(2,sign{j}<0.045),'.');
    result(j) = norm(Re{j}-T,'fro');%norm(Re{j}-project_sphere(Re{j}),'fro');
end


function [steps, data_move, sign] = algorithm(data, data_move, sigma, q, d)
    epsion = 1e-10;
    max_iter = 10000;
    steps = 0;
    n = size(data_move,2);
    sign = zeros(1,n);
    step = 0.5;
    for i = 1:n
 %       for k = 1:max_iter
         k = 0;
         [direction, ~] = mean_shift(data_move(:,i), sigma, data, d, step, q);
         while norm(direction) > epsion
            [direction, sig] = mean_shift(data_move(:,i), sigma, data, d, step, q);
            data_move(:,i) = data_move(:,i)+ direction;
    %        if norm(direction) < epsion
            if k > max_iter
                sign(i) = sig;
                break;
                k
            end  
            k = k+1;
         end
        steps = steps+k/n;
        fprintf('i=%d, sign=%.4f, threshold=%d\n',i, sig, sig<sigma^2/2);
    end
end


function [g,sign_] = mean_shift(x, sigma, data, d, step, q)
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
    BB = B/sum_r-((1-q)*(c-x)*(c-x)');
    [U,E] = eig(BB);
    [~,ind] =sort(diag(E),'descend');
    V = U(:,ind);
    P = eye(size(x, 1))-V(:,1:d)*V(:,1:d)';
    g = step*P*(c - x);
    e = eigs(BB);
    sign_ = e(d+1);
end


function [data_ini, samples, X] = generate_sphere(sigma, D, NumSample, NumIni)
    samples = randn(D, NumSample);
    samples = samples*diag(1./sqrt(sum(samples.^2)))+ sigma*randn(D, NumSample);
    
    data_ini = randn(D, NumIni);
    X = data_ini*diag(1./sqrt(sum(data_ini.^2)));
    data_ini = X + sigma*randn(D,NumIni);
    %0.5*sqrt(sigma)/sqrt(D)*(2*rand(D, NumIni)-1);
end


function [data_ini, samples, X] = generate_random(D, NumSample, NumIni)
    samples = rand(D, NumSample);    
    data_ini = rand(D, NumIni);
    X = data_ini;
end


function re = project_sphere(data)
    re = bsxfun(@rdivide, data, sqrt(sum(data.^2,1)));
end