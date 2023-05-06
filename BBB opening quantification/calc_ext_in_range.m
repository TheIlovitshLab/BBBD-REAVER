function ext_in_segments = ...
    calc_ext_in_range(rcind_seg_cell,...
    all_seg_rads,...
    bw_vessels,...
    redIm,...
    from_px,...
    n_px)
% Function to calculate EB extravasation around deifferent vessel segments
% according to previous vessel segmentation. return statistics for the 
% extravasation in the range (from_px, from_px + n_px] around the vessel
% Inputs:
%     rcind_seg_cell = cell array with every cell containing [row,col]
%         coordinates of a vessel segment
%     all_seg_rads = radii of segments. for dilation.
%     bw_vessels = processed binary image of blood vessles and background
%     redIm = image of EB channel
%     from_px = the starting distance in pixels
%     n_px = how many pixels of extravasation
% Output:
%     eb_ext_in_segments = vector of size n_segments where each element
%        contains the median pixel value of a (j/10)*n_px 
%         neighborhood around the i-th segment in the EB image

ext_in_segments = zeros(length(rcind_seg_cell),1); % Create a placeholder
% Rescale the red image to a range of [0,1]
redIm = rescale(double(redIm));

se_from = strel('disk',from_px,0); 
se_to = strel('disk',n_px,0); 
for n=1:size(rcind_seg_cell,1)  % loop through all segments
    se_small = strel('disk',ceil(all_seg_rads(n)),0); % Create a specific structure element for wire dilation 
    centerline = sub2ind(size(bw_vessels), rcind_seg_cell{n}(:,1),rcind_seg_cell{n}(:,2));    % get all wire-frame elements of the segment (represent the centerline)
    single_seg_bw = false(size(redIm));
    single_seg_bw(centerline) = 1;    % set only the segment wireframe to True
    single_vessel_mask = blockdilate(single_seg_bw,se_small); % Dilate the wire
    single_vessel_mask = single_vessel_mask & bw_vessels; % constrain to only what is inside the vessel
    if from_px == 0
        start_dist = single_vessel_mask;
    else
        start_dist = blockdilate(single_vessel_mask,se_from); % Dilate the single vessel
    end
    end_dist = blockdilate(start_dist, se_to);
    single_vessel_perivasc_mask = logical(end_dist-start_dist);
    single_vessel_perivasc_mask = ...
        logical((single_vessel_perivasc_mask | bw_vessels)-bw_vessels); % Make sure not to take anything thats' within a vessel
    in_vessel_mask = blockdilate(single_seg_bw,strel('disk',ceil(all_seg_rads(n)),0));
    in_vessel_mask = in_vessel_mask & bw_vessels;
    if n_px > 0
        ext_in_segments(n) = ...
            double(median(redIm(single_vessel_perivasc_mask),'all'));
    else
        ext_in_segments(n) = double(median(redIm(in_vessel_mask),'all'));
    end
end
end