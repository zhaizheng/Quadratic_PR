% 
% mu_1 = [-2,-3];
% sigma_1 = 0.3*[0.3 1; 1 3];
% data_1 = mvnrnd(mu, sigma, 500);
% 
% mu_2 = [3,2];
% sigma_2 = 0.3*[4 1.5; 1.5 1];
% data_2 = mvnrnd(mu_2, sigma_2, 500);
% 
% plot(data_1(:,1), data_1(:,2),'+')
% hold on
% plot(data_2(:,1), data_2(:,2),'o')
function data = data_generation(step,var)
    t = 3:step:9;
    x = t.^(0.4);
    y = var.*randn(size(t));
    %y = sin(t)+var.*sqrt(t/3).*randn(size(t));
    x2 = t;
    y2 = 2*cos(t)+t+var*randn(size(t))+3;
%     figure
%     plot(x,y,'*')
%     hold on
%     plot(x2, y2,'ro')
%     X = [x, x2];
%     Y = [y, y2];
    data=[x;y];
end

    