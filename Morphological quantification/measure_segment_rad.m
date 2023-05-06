function all_seg_rads =...
    measure_segment_rad(rcind_seg_cell,bw_seg)
%{
    Function to measure the mean radius of each segment
    Input arguements:
        rcind_seg_cell = vertical cell array. every cell contains the 
                            [row,col] ccordinates of all the wire-frame
                            pixels belonging to that segment
        bw_seg = BW image of segmented vessels
        rc_ep = [row,col] coordinates of all end points
    Output arguements:
        all_seg_rads = table with number of rows equal to number of
            segments, every row contains the median and max redius of
            that segment
%}
    
d_nn = bwdist(~bw_seg);    % for each point calculate the distance to the nearest black point (background).
   
max_seg_rads = zeros(1, size(rcind_seg_cell,1));
median_seg_rads = zeros(1, size(rcind_seg_cell,1));

for n=1:size(rcind_seg_cell,1)  % loop through all segments
    lind_seg = sub2ind(size(bw_seg), rcind_seg_cell{n}(:,1),rcind_seg_cell{n}(:,2));    % get all wire-frame elements of the segment (represent the middle-line)
    max_seg_rads(n) = max(d_nn(lind_seg));    
    median_seg_rads(n) = median(d_nn(lind_seg));
    % Measure circularity
    
end
all_seg_rads = ...
    table(max_seg_rads',median_seg_rads',...
    'VariableNames',{'max','median'});
end

