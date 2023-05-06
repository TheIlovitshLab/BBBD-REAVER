function results_tbl = batch_analyze_vessel_density(path)
% Function to iterate over a directory of images pre-processed by REAVER, with 
% Existing ".mat" analysis files, and compute vessel density
% Input:
%     path = path to a folder pre-processed by REAVER
% Output arguements:
%     results_tbl = table of n rows corresponding to n files in path. for
%         each file it contains the vessel density

    %% Predefined params
    if nargin < 1
        path = uigetdir();
    end
    %% Load files
    all_files = dir(fullfile(path,'*.mat'));
    summary_idx = cellfun(@(x) strcmp(x,'User Verified Table.mat'),...
        {all_files.name});
    summary = load(fullfile(path,all_files(summary_idx).name));
    data_files = all_files(~summary_idx);
    EB_analysis = cellfun(@(x) strfind(x,'EB_analysis'),{data_files.name},'UniformOutput',0);
    non_analysis = cellfun(@(x) isempty(x),EB_analysis);
    data_files = data_files(non_analysis);
    verified = summary.userVerified(cell2mat(summary.userVerified(:,2)),1);
    verified = cellfun(@(x) strrep(x,'.tif','.mat'),verified,'UniformOutput',false);
    [~,~,verified_idx] = intersect(verified',{data_files.name});
    verified_files = data_files(verified_idx);   % Get only verified files

    %% create placeholders
    n_files = numel(verified_files);
    image_name = {verified_files.name}';
    vessel_density = cell(n_files,1);
    length_density = cell(n_files,1);
    bifurcation_density = cell(n_files,1);
    results_tbl = ...
        table(image_name,vessel_density, length_density,bifurcation_density);
    %% Iterate over files
    for i = 1:n_files
        fprintf('Processing file %d of %d\n',i,numel(verified_files));
        % Calc extravasation for every segment.
        metric_st.vessel_density = cell(1,1);
        metric_st.vessel_density{:} = ...
            v_density(fullfile(path,verified_files(i).name));
        metric_st.length_density = cell(1,1);
        metric_st.length_density{:} = ...
            l_density(fullfile(path,verified_files(i).name));
        metric_st.bifurcation_density = cell(1,1);
        metric_st.bifurcation_density{:} = ...
            b_density(fullfile(path,verified_files(i).name));
        metric_st.image_name = image_name(i);
        results_tbl(i,:) = struct2table(orderfields(metric_st,[4,1,2,3]));
    end
end