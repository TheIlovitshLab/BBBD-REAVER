function combine_image_with_edges(multi_line, im_path)
% Function to write a new 2-channel tiff where the segmentation outline is
% overlayed on the original image
% Inputs:
%   im_path = (optional) path to the tiff file
%   multi_line = (optional) array of distances in pixels to also draw
%       around each segmented vessel
% Output:
%   saves a copy of the original image with the overlayed otline at the
%   same file location and the same file name with a "_combined" suffix
if nargin < 2
    [im_file, folder] = uigetfile('*.tif');
    im_path = fullfile(folder,im_file);
end
redIm = imread(im_path,1);
greenIm = imread(im_path,2);
composite = cat(3,redIm,greenIm);
matfile = load(strrep(im_path,'.tif','.mat'));
original_bw = matfile.derivedPic.BW_2;
eroded = imerode(original_bw,strel('disk',1));
edges = original_bw-eroded;
edges_3d = logical(repmat(edges,1,1,2));
composite(edges_3d) = 1023;
if nargin > 0  % User specified multiline
    for i = 1:length(multi_line)
        dilated = imdilate(original_bw,strel('disk',multi_line(i)));
        eroded = imerode(dilated,strel('disk',2));
        edges = dilated-eroded;
        edges_3d = logical(repmat(edges,1,1,2));
        composite(edges_3d) = 1023;
    end
end
rgb_comp = 64*composite;
rgb_comp(:,:,3) = 0;
imagesc(rgb_comp);
bfsave(composite,strrep(im_path,'.tif','_combined.tiff'));
end