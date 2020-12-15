function [X, X1] = sphere2(sigma,num)
    data = mvnrnd([0,0,0], eye(3), num)';
    X = bsxfun(@rdivide, data, sqrt(sum(data.^2,1)));
    %X1 = bsxfun(@times, X, 1+sigma*randn([1,num]));
    X1 = X + sigma*randn([3,num]);
end