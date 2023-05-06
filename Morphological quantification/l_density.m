function vessel_density = l_density(mat_path)
% Function to extract the vessel length density [1/mm]
% Input:
%     mat_math = path to a .mat file containing the verified image
%                parameters, BW image and wireframe
st = load(mat_path);
image_area = st.imageSize(1)*st.imageSize(2)*(st.image_resolution^2);
vessel_density = (st.metrics.vesselLength*st.image_resolution)/image_area;
vessel_density = vessel_density*10^3; % To go from um^-1 to mm^-1
end