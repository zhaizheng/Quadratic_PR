nsample = 4000;
data = randn(2, nsample);
ndata = normal(data);

ndata3 = [ndata;zeros(1,nsample)]+0.03.*randn(3,nsample);
cdata3 = [normal(ndata3(1:2,:)); zeros(1,nsample)];
fdata = 0.2*normal(ndata3 - cdata3)+0.03.*randn(3,nsample) +cdata3;
scatter3(fdata(1,:),fdata(2,:),fdata(3,:))

function ndata = normal(data)
    ndata = data*diag(1./sqrt(sum(data.^2, 1)));
end