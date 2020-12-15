[images,labels] = read_file();
interest = images(:,:,labels==8);

function display(interest)
    num = 7;
    I = zeros(28*num,28*num);
    for i = 1:num
        for j = 1:num
            I((i-1)*28+1:i*28, (j-1)*28+1:j*28) = interest(:,:,(i-1)*10+j);
        end
    end
    imshow(I)
end