

%imshow(reshape(Data(:,1),64,[]));

%img = (rem(randperm(64*64),256))';
obj = {'obj1_*','obj13*'};
for k = 1:2
    WHOLE = [];
    IND = [];
    Data = [];
    dd = fullfile(obj{k});
    dirOutput=dir(fullfile(obj{k}));
    name = dirOutput.name;

    for i = 1:length(dirOutput)
        temp = imread(dirOutput(i).name);
        vector = reshape(temp(1:2:end,1:2:end),[],1);
        Data = [Data, vector];
    end
    for rep = 1:5
        img = floor(256*rand(64*64,1));
        for step = 1:10000
            [direction,ind] = Hc77(img, 1e4, double(Data), 1, 0.3);
            if norm(direction)> 0.1
                img = img + direction;
            else
                break;
            end
            fprintf('step=%d,norm=%f\n',step,norm(direction));
        end
        S_Data = Data(:,ind(1:5));
        TEMP = [];
        for i = 1:5
            TEMP = [TEMP,reshape(Data(:,ind(i))',64,[])];
        end
        IND = [IND;ind];
        LINE = [TEMP,uint8(reshape(img',64,[]))];
        WHOLE = [WHOLE; LINE];
    end
    subplot(1,2,k);
    imshow(WHOLE);
end



function [g,ind] = Hc77(x, sigma, data, d, step)

    sq_distance = sum((data - x).^2,1);
    [~, ind] = sort(sq_distance);
    
    c = zeros(size(x));
    sum_r = 0;
    %for i = 1:size(data,2)
    s = 5;
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
    [V,D] = eigs(B);
    P = eye(size(x, 1))-V(:,1:d)*V(:,1:d)';
    %g = 2*sum_r/size(data,2)/sigma^2*P*(c - x);
    g = step*P*(c - x);
end


