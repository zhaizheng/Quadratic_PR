x = linspace(-10,10,100);
y = sin(x)+0.1*randn(1,100);

samples = [x;y];
plot(x, y,'*');
hold on
test = [(max(x)-min(x))*(rand(1,1000)-0.5+(max(x)+min(x))/2);(max(y)-min(y))*(rand(1,1000)-0.5+(max(y)+min(y))/2)];
sigma = 1;
for i =1:500
    c1 = C(samples, sigma, test(:,i));
    c2 = C(samples, sigma, c1);
    c3 = C(samples, sigma, c2);
    c4 = C(samples, sigma, c3);
    plot([c3(1),c4(1)],[c3(2),c4(2)],'-')
    hold on
end


function c = C(data, sigma, x)
    c = zeros(size(x));
    sum_r = 0;
    for i = 1:size(data,2)
        r = exp(-norm(x - data(:,i)).^2/(sigma^2));
        c = c + r*data(:,i);
        sum_r = sum_r + r;
    end
    c = c/sum_r;
end