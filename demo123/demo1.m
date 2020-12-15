t = -1:0.2:1;
y = sqrt(1- 1/3*t.^2)+0.06*randn(size(t));


X = [t;y];

x0a = [1;0.25];

for i = 1:3
    subplot(1,3,i); 
    subplot('position', [0.04*i+0.29*(i-1) 0.06 0.29 0.9]);
    plot(t,sqrt(1- 1/3*t.^2),'--','LineWidth',3);
    hold on
    plot(t, y, '*','MarkerSize',10);
    hold on
    x0 = [0; i*0.3];
    generate_figure(X, x0, x0a);
    set(gca,'FontSize',14);
end


function generate_figure(X, x0, x0a)

    plot(x0(1), x0(2), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'r');

    [U,L] = principal(X, x0);

    hold on
    c1 = positive(U(:,1),x0a);
    c2 = positive(U(:,2),x0a);
    h1 = quiver(x0(1),x0(2),c1*L(1)*U(1,1),c1*L(1)*U(2,1),'linewidth',2,'color',[1,0,0]);
    set(h1,'AutoScale','off', 'AutoScaleFactor', 2);
    hold on
    h2 = quiver(x0(1),x0(2),c2*L(2)*U(1,2),c2*L(2)*U(2,2),'linewidth',2,'color',[0,0,1]);
    set(h2,'AutoScale','off', 'AutoScaleFactor', 2);
%     hold on
%     t = 0:0.1:pi;
%     x = 0.8*cos(t);
%     y = 0.8*sin(t)+0.3;
% 
%     plot(x,y,'-','LineWidth',2);
    axis([-1, 1.1, 0.1, 1.3]);
end


function c = positive(a,b)
    if a'*b>0
        c = 1;
    else
        c = -1;
    end
end


function [U2,l2] = principal(X, x0)
    S = zeros(size(X, 1));
    for i=1:size(X,2)
        S = S + exp(-norm(X(:,i)-x0)^2/0.4)*(X(:,i)-x0)*(X(:,i)-x0)';
  %      S = S + (X(:,i)-x0)*(X(:,i)-x0)';
    end
%    [U2,L] = svd(S);
%    l2 = diag(L);
    [U,L] = eig(S);
    [l2,ind] = sort(diag(L),'descend');
    %l2 = l(ind);
    U2 = U(:,ind);
end