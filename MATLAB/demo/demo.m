
[real, samples] = generate_data();

subplot(1,2,1)

plot(samples(1,:),samples(2,:),'*');
hold on
plot(real(1,:),real(2,:),'-');
sigma = 1;
for i = 1:size(samples,2)
    e = Hessian(real(:,i),samples, sigma);
    quiver(real(1,i),real(2,i),e(1),e(2));
    hold on
end

subplot(1,2,2)

plot(samples(1,:),samples(2,:),'*');
hold on
plot(real(1,:),real(2,:),'-');

for i = 1:size(samples,2)
    e = Hessian2(real(:,i),samples, sigma);
    quiver(real(1,i),real(2,i),e(1),e(2));
    hold on
end


function [real, samples] = generate_data()
    theta = -2*pi:0.1:2*pi;
    x = theta;
    n = 0.08*randn(1,length(x));
    y = sin(theta)+n;
    z = sin(theta);
    samples = [x;y];
    real = [x;z];
%     plot(x,y,'*');
%     hold on
%     plot(x,z,'-');
end

function e = Hessian2(x, samples, sigma)
    H = zeros(length(x));
    for i = 1:size(samples,2)
        H = H + exp(-norm(x-samples(:,i)).^2/(sigma^2))*(x-samples(:,i))*(x-samples(:,i))';
    end
    [U,~,~] = svd(H);
    e = U(:,1);
end

function e = Hessian(x, samples, sigma)
    H = zeros(length(x));
    sq_distance = sum((samples - x).^2,1);
    [~, ind] = sort(sq_distance);
    for i = ind(1:15)%size(samples,2)
        H = H + exp(-norm(x-samples(:,i)).^2/(sigma^2))*(x-samples(:,i))*(x-samples(:,i))';
    end
    [U,~,~] = svd(H);
    e = U(:,1);
end