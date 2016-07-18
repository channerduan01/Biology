function displayImages(dataset, original_image_size, zoom_scale, str_title)
    pad = 1;
    real_image_size = floor(original_image_size/zoom_scale);
    dim = size(dataset,2);
    display_rows = floor(sqrt(dim));
    display_cols = ceil(dim / display_rows);
    display_array = -1 * ones(pad + display_rows * (real_image_size(1) + pad), pad + display_cols * (real_image_size(2) + pad));
    for j = 1:display_rows
        for i = 1:display_cols
            index = (j-1)*display_cols+i;
            if index > dim, break; end;
            comp = imresize(reshape(dataset(:,index), original_image_size), real_image_size);
            comp = normalize(comp);
            display_array(pad + (j - 1) * (real_image_size(1) + pad) + (1:real_image_size(1)), ...
                pad + (i - 1) * (real_image_size(2) + pad) + (1:real_image_size(2))) = comp;   
        end
    end
    figure();
    colormap(gray);
    imagesc(display_array,[-1,1]);
    axis('off');
    axis('equal');
    title(str_title, 'FontSize', 20);
end