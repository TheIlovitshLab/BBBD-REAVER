function bifurcation_density = b_density(mat_path)
% Function to extract the number of bifurcations per area
%     mat_math = path to a .mat file containing the verified image
%                parameters, BW image and wireframe
st = load(mat_path);
n_bif = length(st.derivedPic.branchpoints);
image_area = st.imageSize(1)*st.imageSize(2)*(st.image_resolution^2);
bifurcation_density = n_bif/image_area;
bifurcation_density = bifurcation_density*10^6; % To go from um^-2 to mm^-2
end