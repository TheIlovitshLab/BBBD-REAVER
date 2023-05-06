density_tbl = batch_analyze_vessel_density;
im_names_no_suffix = cellfun(@(x) x(1:end-4), density_tbl.image_name,...
    'UniformOutput',false);
regions =cellfun(@(x)region_from_name(x), im_names_no_suffix,...
    'UniformOutput',false);
unique_regions =  unique(regions);
vessel_density_per_region = cell(numel(unique_regions),1);
length_density_per_region = cell(numel(unique_regions),1);
bifurcation_density_per_region = cell(numel(unique_regions),1);

for i = 1:numel(unique_regions)
   vessel_density_per_region{i} = cell2mat(density_tbl.vessel_density(...
       strcmp(unique_regions{i},regions)~=0));
   length_density_per_region{i} = 3*cell2mat(density_tbl.length_density(...
       strcmp(unique_regions{i},regions)~=0));
   bifurcation_density_per_region{i} = 3*cell2mat(density_tbl.bifurcation_density(...
       strcmp(unique_regions{i},regions)~=0));
end
new_tbl = table(unique_regions,...
    vessel_density_per_region,...
    length_density_per_region,...
    bifurcation_density_per_region);
%% Plotting
% for i = 1:numel(unique_regions)
%     Violin(new_tbl.bifurcation_density_per_region(i),i,...
%         'ViolinColor',{[1,0,0]});
% end
boxplot2(new_tbl.bifurcation_density_per_region)
xticks(1:numel(unique_regions));
xticklabels(unique_regions);
%%
v_tbl = new_tbl(:,[1,2]);
l_tbl = new_tbl(:,[1,3]);
b_tbl = new_tbl(:,[1,4]);
%% Saving
folder = uigetdir();
group = 'test';  % Change accordingly
writetable(v_tbl,fullfile(folder,[group,'_v.csv']));
writetable(l_tbl,fullfile(folder,[group,'_l.csv']));
writetable(b_tbl,fullfile(folder,[group,'_b.csv']));
%% Helper function
function region = region_from_name(name)
% Helper function to retrieve the brain region (3rd arguement) from an
% image name with the regular expression
% "<animal>_<section>_<region>_<index>
sp = textscan(name,'%s','Delimiter','_');
region = sp{1}{3};
end