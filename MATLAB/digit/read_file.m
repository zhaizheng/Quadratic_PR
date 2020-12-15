function [images,labels] = read_file()
    images = load_images();
    labels = load_labels();
end
function images = load_images()
    filename = 'train-images-idx3-ubyte';
    FID = fopen(filename,'r');
    magic = fread(FID, 1 , 'int32', 0, 'ieee-be');
    numImages = fread(FID, 1 , 'int32', 0, 'ieee-be');
    numRows = fread(FID, 1 , 'int32', 0, 'ieee-be');
    numCols = fread(FID, 1 , 'int32', 0, 'ieee-be');
    images = fread(FID, inf, 'unsigned char');
    images = reshape(images, numCols, numRows, numImages);
    images = permute(images, [2,1,3]);
    fclose(FID)
end


function labels = load_labels()
    filename = 'train-labels-idx1-ubyte';
    fid = fopen(filename, 'r');
    magic = fread(fid, 1, 'int32', 0, 'ieee-be');
    numLabels = fread(fid, 1, 'int32', 0, 'ieee-be');
    labels = fread(fid, Inf, 'unsigned char');
    labels = labels';
    fclose(fid);
end
