
[real, samples] = generate_data();

subplot('position',[0.04 0.06 0.45 0.9])

plot(samples(1,:),samples(2,:),'*');
hold on
plot(real(1,:),real(2,:),'-','linewidth',2);
sigma = 1.5;

num = floor(size(samples,2)/2);

for i = 2*(1:num)
    e = Hessian(real(:,i),samples, sigma);
    h0 = quiver(real(1,i),real(2,i),e(1),e(2),'linewidth',2);
    set(h0, 'AutoScale','off','AutoScaleFactor',2)
    hold on
end
axis([-7.4,7,-1.5, 1.5])
set(gca,'FontSize',14);
subplot('position',[0.53 0.06 0.45 0.9])

plot(samples(1,:),samples(2,:),'*');
hold on
plot(real(1,:),real(2,:),'-');

%num = floor(size(samples,2)/3);
for i = 2*(1:num)
    [e,l] = Hessian2(real(:,i),samples, sigma);
    h = quiver(real(1,i),real(2,i),e(1),e(2),'linewidth',2);
    set(h, 'AutoScale','off','AutoScaleFactor',2)
    hold on
end
axis([-7.4,7,-1.5, 1.5])
set(gca,'FontSize',14);

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

function [e,l] = Hessian2(x, samples, sigma)
    H = zeros(length(x));
    for i = 1:size(samples,2)
        H = H + exp(-norm(x-samples(:,i)).^2/(sigma^2))*(x-samples(:,i))*(x-samples(:,i))';
    end
    [U,L,~] = svd(H);
    e = U(:,1);
    l = L(1,1); 
end

function e = Hessian(x, samples, sigma)
    H = zeros(length(x));
    sq_distance = sum((samples - x).^2,1);
    [~, ind] = sort(sq_distance);
    for i = ind(1:20)%size(samples,2)
        H = H + exp(-norm(x-samples(:,i)).^2/(sigma^2))*(x-samples(:,i))*(x-samples(:,i))';
    end
    [U,~,~] = svd(H);
    e = U(:,1);
end