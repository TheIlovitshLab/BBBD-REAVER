function metric_st = features_from_frame(mat_path,n_px)
% Custom image processing and measurements extraction function.
% Inputs:
%     mat_math = path to a .mat file containing the verified image
%                parameters, BW image and wireframe
%     n_px = Size of neighborhood around a blood vessel (in px)
% Output:
%     metric_st = structure with all measurements. containing the following
%         fields:                   
%             'vessel_length_um'            
%             'segments_diam_um'
%             'bbbd_marker_median'

                   
warning('off','imageio:tiffmexutils:libtiffWarning');

if isempty(dir(mat_path))
    error('File not found in specified path: %s\n' , mat_path);
end

% Load matlab reaver file
st = load(mat_path);
umppix = st.image_resolution;

% Load red channel
im_file = [mat_path(1:end-4),'.tif'];
t = Tiff(im_file);
setDirectory(t,1)
redIm = read(t);   % We only need the red channel
% Initialize output struct
metric_st=struct();

% Add average radius and calssify each lineseg
rcind_seg_cell = skel_2_linesegs(st.derivedPic.wire,...
    fliplr(st.derivedPic.branchpoints),fliplr(st.derivedPic.endpoints));

% Measure segment radii and record diameter
all_seg_rads = measure_segment_rad(rcind_seg_cell, st.derivedPic.BW_2);
%Multiply by 2 and add the pixel for the center point to get the diameter
median_diams = 2.*all_seg_rads.median+1;
max_diams = 2.*all_seg_rads.max+1;

median_segment_diam_um = cell(1,1);
max_segment_diam_um = cell(1,1);

median_segment_diam_um{1,1} = median_diams .* umppix;
max_segment_diam_um{1,1} = max_diams .* umppix;

metric_st.median_segment_diam_um = median_segment_diam_um;
metric_st.max_segment_diam_um = max_segment_diam_um;

% Measure segments lengths in um
all_segment_len_um = cell(1,1);
all_segment_len_um{1,1} = cellfun(@(x) size(x,1),rcind_seg_cell).*umppix;
metric_st.segment_len_um = all_segment_len_um;

bbbd_marker_median = cell(1,1);
bbbd_marker_median{1,1} = ...
    calc_ext_in_range(rcind_seg_cell,...
    all_seg_rads.max,st.derivedPic.BW_2,...
    redIm,...
    n_px);

metric_st.bbbd_marker_median = bbbd_marker_median;

% Remove nans
nans = isnan(bbbd_marker_median{1});
metric_st.bbbd_marker_median{1}(nans) = [];
metric_st.median_segment_diam_um{1}(nans) = [];
metric_st.max_segment_diam_um{1}(nans) = [];
metric_st.segment_len_um{1}(nans) = [];
end

