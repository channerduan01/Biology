function [images, subject_num, samples_num, size_img, lables] = dataLoad(database,isFlat)
    subject_num = length(database);
    samples_num = database(1).Count;
    rows = size(read(database(1),1),1);
    columns = size(read(database(1),1),2);
    size_img = [rows, columns];
    if isFlat
        images = zeros(subject_num*samples_num, rows*columns);
    else
        images = zeros(rows, columns, subject_num*samples_num);
    end
    lables = zeros([1, subject_num*samples_num],'int8');
    index = 1;
    for i=1:subject_num
        for j=1:samples_num
            if isFlat
                images(index, :) = reshape(read(database(i),j), 1, []);
            else
                images(:, :, index) = read(database(i),j);
            end
            lables(index) = i;
            index = index+1;
        end
    end
end