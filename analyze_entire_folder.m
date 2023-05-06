function file_path = ...
    analyze_entire_folder(n_px, path, from_px)
% Function to iterate over a directory of images pre-processed by REAVER 
% (with pre-existing ".mat" analysis files), and compute BBBD marker leakage
% Input arguements:
%     n_px = width of perivascular are in pixels
%     path = path to a folder pre-processed by REAVER
%     from_px = distance from vessel wall to start the eb extraction from
%       [pixels]
% Output arguements:
%     file_path - path to a .mat file that will contain a table of n rows
%         corresponding to n files in path. for each file it contains all
%         vessel segments median diameter and extravasation (median BBBD 
%         marker intensity in perivascular area). The filename also 
%         contains the processing hyperparameters from_px and n_px

    %% Predefined params
    if nargin < 3
        from_px = 0;
    end
    if nargin < 2
        path = uigetdir();
    end
    if nargin == 0
        n_px = 10;
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
    median_segment_diam_um = cell(n_files,1);
    max_segment_diam_um = cell(n_files,1);
    bbbd_marker_median = cell(n_files,1);
    segment_len_um = cell(n_files,1);
    results_tbl = ...
        table(image_name,median_segment_diam_um,...
        max_segment_diam_um,segment_len_um,...
        bbbd_marker_median);
    %% Iterate over files
    parfor i = 1:n_files
        fprintf('Processing file %d of %d\n',i,numel(verified_files));
        % Calc extravasation for every segment.
        metric_st = ...
            features_from_frame(fullfile(path,verified_files(i).name),...
            n_px, from_px);
        metric_st.image_name = image_name(i);
        results_tbl(i,:) = ...
            struct2table(orderfields(metric_st,[5,1,2,3,4]));
    end
    res = struct('table',results_tbl,'n_px',n_px,'from_px',from_px);
    str = 'Ext_analysis_';
    if from_px ~= 0
        str = [str,'from_',num2str(from_px),'px_'];
    end
    str = [str,'_to_',num2str(from_px+n_px),'px'];
    str = [str,'.mat'];
    file_path = fullfile(path,str);
    save(file_path,'res');
end