% This script requires manual loading of a .mat file of a segmented image
%% If you want a bright background and black segmentation
empty_mask = false(1200,1920);
comp = 255*repmat(~derivedPic.BW_2,1,1,3);
wire = imdilate(derivedPic.wire,strel('disk',2,0));
% If you want a blue wire replace the vessel mask to the 3rd channel
vessel_mask = repmat(empty_mask,1,1,3);
vessel_mask(:,:,2) = wire;
comp(vessel_mask) = 255;
%% If you want a black background and white segmentation
empty_mask = false(1200,1920);
comp = 255*repmat(derivedPic.BW_2,1,1,3);
wire = imdilate(derivedPic.wire,strel('disk',2,0));
% If you want a blue wire replace the vessel mask to the 3rd channel
vessel_mask = repmat(wire,1,1,3);
vessel_mask(:,:,2) = empty_mask;
comp(vessel_mask) = 0;

%% Plotting
imshow(comp)
%% Saving
imwrite(comp,...
    fullfile("C:\Users\Admin\Documents\Studies\Masters\Research\My papers\Sharon's paper folder\Figures\image processing scheme",...
    'segmentation and wire in green.jpg'));