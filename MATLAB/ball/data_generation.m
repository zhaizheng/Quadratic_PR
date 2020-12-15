function [X, X1] = sphere()
    alpha = 0:0.3:pi;
    theta = 0:0.3:pi;

    [Alpha, Theta] = meshgrid(alpha,theta);
    epsilon = 0.1*randn(size(Alpha));


    x = cos(Alpha).*sin(Theta);
    y = cos(Alpha).*cos(Theta);
    z = sin(Alpha);

    x1 = (1.+epsilon).*cos(Alpha).*sin(Theta);
    y1 = (1.+epsilon).*cos(Alpha).*cos(Theta);
    z1 = (1.+epsilon).*sin(Alpha);

    % [x2,y2,z2]=sphere(1000);

    X = [x; y; z];
    X1 = [x1; y1; z1];

% scatter3(x(:),y(:),z(:), '*')
% hold on
% mesh(x1,y1,z1)
end


