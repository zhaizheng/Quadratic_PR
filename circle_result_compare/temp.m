figure
t = tiledlayout(2,3,'TileSpacing','Compact');
nexttile
plot(X2(1,1:1000),X2(2,1:1000),'.','MarkerSize',10);
hT = title('Original Data','interpreter','tex');
set(hT, 'FontSize', 14)
%Title = {'log(x)','p^{1/4}(x)','p^{1/2}(x)','p^{3/4}(x)','p(x)'};
Title = {'\Gamma(-10000,x)','\Gamma(-1,x)','\Gamma(0,x)','\Gamma(0.5,x)','\Gamma(1,x)'};
for j = 1:5
    nexttile
    plot(T(1,:),T(2,:),'.');
    hold on
    plot(Re{j}(1,sign{j}<0.045),Re{j}(2,sign{j}<0.045),'.','MarkerSize',10);
    hT = title(Title{j},'interpreter','tex');
    axis([-2, 2, -2 2]);
    %legend({'truth', 'ridge'})
    %set(gca,'fontsize',14)
    set(hT, 'FontSize', 14)
end
