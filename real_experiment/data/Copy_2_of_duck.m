%imshow(reshape(Data(:,1),64,[]));

%img = (rem(randperm(64*64),256))';
d = 50;
file_name = {'obj1_*','obj13*'};


WHOLE = [];
IND = [];
Data = extract_data(file_name{1});
[lD, Us] = low_dimension(Data, d);
sigma = 0.2*norm(lD,'fro')/size(Data,2);
noisy_lD = lD + sigma*randn(size(lD));

error = zeros(4,5,4);
output = cell(4,5,4);
pic_idx = 1:10;

Neig = 14:2:20;
for i = 1:length(Neig)
    for rep = 1:length(pic_idx)
    %    img = floor(256*rand(d,1));
        img = noisy_lD(:,pic_idx(rep));
        data = noisy_lD; step = 0.3; d = 1; sigma = 1e4; beta = 2; %neig = 4;

        for alg = 1:5
            output{alg,rep,i}= run_algorithm(img, data, d, step, alg, Neig(i));
            error(alg,rep,i) = norm(lD(:,rep)-output{alg, rep, i},'fro');
        end

        %new_imag = [ld_to_img(Us,lD(:,rep)),ld_to_img(Us,noisy_lD(:,rep)), ld_to_img(Us,img)];
        %[uint8(reshape((Us*img)',128,[])),uint8(reshape((Us*img)',128,[])),uint8(reshape((Us*lD(:,rep))',128,[]))];
        %WHOLE = [WHOLE; new_imag];
        %WHOLE = [WHOLE; reconstruct(img, Us, noisy_lD, ind)];
    end
end
%%
%subplot(1,6,1)
%subplot('Position',[0.03 0.1 0.1,0.85])
figure(1)
t = tiledlayout(1,5,'TileSpacing','Compact');
nexttile
Orig = [];
for rep = 1:6
    row = [ld_to_img(Us,lD(:,rep)),ld_to_img(Us,noisy_lD(:,rep))];
    Orig = [Orig; row];
end
imshow(Orig)

%Ti = {'{\it l}-SCRE','SCRE', 'MFIT-ii','MFIT-i'};

for neig = 1:4
    nexttile
    imag = [];
    for rep = 1:6
        row = [];
        for alg = 1:5
            row = [row,ld_to_img(Us,output{alg,rep,neig})];
        end
        imag = [imag; row];
    end
    imshow(imag);
    title(['Neig=',num2str(Neig(neig))]);
end
   
figure(2)
%subplot('Position',[(0.03+0.15)*5-0.04 0.1 0.12,0.85])
marker = {'-o','-^','-*','-d','->'};
for alg = 1:5
    plot(Neig, mean(squeeze(error(alg,:,:)),1),marker{alg},'linewidth',2);
    hold on
end 
legend('{\it l}-SCRE:q=1','{\it l}-SCRE:q=0.8', '{\it l}-SCRE:q=0.6','MFIT-ii','MFIT-i');
set(gca,'FontSize',12);
%%
function obj = run_algorithm(obj, data, d, step, alg, k)
    dis_sort = sort(sqrt(sum((data-obj).^2,1)),'ascend');
    sigma = dis_sort(k); r = dis_sort(k);
    beta = 2;
    neig = 18; 
    for rep = 1:10000
        switch alg
            case 1
                q = 1;
                direction = local_mean_shift_parametrized(obj, sigma, data, d, step, neig, q);
            case 2
                q = 0.8;
                direction = local_mean_shift_parametrized(obj, sigma, data, d, step, neig, q);
            case 3
                q = 0.6;
                direction = local_mean_shift_parametrized(obj, sigma, data, d, step, neig, q);
            case 4
                direction = xia(obj, data, r, beta, d, step);
            case 5
                direction = mfit(obj, data, r, beta, d, step);
        end
        if norm(direction)> 0.01
            obj = obj + direction;
        else
            break;
        end
        %fprintf('step=%d,norm=%f\n',step,norm(direction));
    end
end

function Img = ld_to_img(Us, vec)
    Img = uint8(reshape((Us*vec)',128,[]));
end


function LINE = reconstruct(img, Us, Data, ind)
        TEMP = [];
        for i = 1:5
            TEMP = [TEMP,reshape(Data(:,ind(i))',128,[])];
        end
        LINE = [TEMP,uint8(reshape((Us*img)',128,[]))];
end

function [lD, Us] = low_dimension(Data, d)
    Data = double(Data);
    [U,~,~] = svd(Data);
    lD = U(:,1:d)'*Data;
    Us = U(:,1:d);
end



function Data = extract_data(name)
    Data = [];
    dirOutput=dir(fullfile(name));
    name = dirOutput.name;

    for i = 1:length(dirOutput)
        temp = imread(dirOutput(i).name);
        vector = reshape(temp,[],1);
        Data = [Data, vector];
    end
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



function g = Hc77(x, sigma, data, d, step)

    sq_distance = sum((data - x).^2,1);
    [~, ind] = sort(sq_distance);
    
    c = zeros(size(x));
    sum_r = 0;
    %for i = 1:size(data,2)
    s = 15;
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
    [V,~] = eigs(B);
    P = eye(size(x, 1))-V(:,1:d)*V(:,1:d)';
    %g = 2*sum_r/size(data,2)/sigma^2*P*(c - x);
    g = step*P*(c - x);
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

