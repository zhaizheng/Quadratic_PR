%%
marker = {'o','>','d','*'};
% scatter(outputs{1}(1,:), outputs{1}(2,:),'bo')
%%
X = randn(2,500);
Xn = X*diag(1./sqrt(sum(X.^2, 1)));
scatter(Xn(1,:), Xn(2,:),'.')
hold on
%%
for i = 7:9
    hold on
    scatter(outputs{i}(1,:), outputs{i}(2,:), marker{i-6})
end