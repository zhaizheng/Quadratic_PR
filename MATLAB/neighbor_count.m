function neighbor_count()
    data = data_generation(0.02,0.3);
    size_ = size(data,2);


    XL = -8:0.1:3;
    YL = -4:0.1:6;
    [SX,SY] = meshgrid(XL,YL);

    Z = zeros(size(SX));
    for i=1:size(SX,1)
        for j = 1:size(SX,2)
            for k = 1:size_
                Z(i,j) = Z(i,j) + exp(-norm([SX(i,j);SY(i,j)]-data(:,k)).^2/0.05) ;
            end
        end
    end
    figure
    imagesc(Z)


    C = zeros(length(XL)*length(YL),size_);
    for i = 1:length(XL)*length(YL)
        for j = 1:size_
            %if norm(data(:,i)-data(:,j)) < 0.01
            C(i,j) = exp(-norm([SX(i);SY(i)]-data(:,j)).^2) ;
            %end
        end 
    end
    neighbor_size = sum(C,2);

    M = zeros(length(XL)*length(YL),15);

    for i = 1:length(XL)*length(YL)
        M(i,:) = base_Calulate(SX(i),SY(i));
    end

    [U,L,V] = svd(M);

    alpha = V'*diag(1./diag(L))*U(:,1:15)'*neighbor_size;



    Z = zeros(size(SX));
    for i = 1:length(XL)
        for j = 1:length(YL)
            Z(j,i) = base_Calulate(SX(j,i),SY(j,i))*alpha;
        end
    end
    figure
    imagesc(Z)
    % hold on
    % plot(data(1,:),data(2,:),'o')
end

