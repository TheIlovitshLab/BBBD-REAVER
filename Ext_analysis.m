classdef Ext_analysis
    % Extravasation analaysis results
    % includes both control and MB + FUS data
    properties
        segment_tbl
        n_px
        from_px
        UM_PX    
    end
    methods
        %% General tabular functions
        function obj = Ext_analysis(varargin)
            % Input arguements:
            %   control_file = control "Ext_analysis_entire_folder" file
            %   test_file = test "Ext_analysis_entire_folder" file
            %   um_px = ratio of microns to pixel (default = 0.29288)
            P = inputParser();
            P.addOptional('control_file',[],@(x) isfile(x));
            P.addOptional('test_file',[],@(x) isfile(x));
            P.addOptional('um_px',0.29288,@(x) isnumeric(x));
            P.parse(varargin{:});
            control_file = P.Results.control_file;
            test_file = P.Results.test_file;
            obj.UM_PX = P.Results.um_px;
            if isempty(control_file)
                % Open file dialogue to select the analysis file of control
                % brains
                [file1,folder1] = uigetfile('*.mat','Choose control analysis file');
                control_file = fullfile(folder1,file1);
            end
            control = load(control_file);
            if isempty(test_file)
                % Open file dialogue to select the analysis file of MB+FUS
                % treated brains
                [file2,folder2] = uigetfile('*.mat','Choose treatment analysis file');
                test_file = fullfile(folder2,file2);
            end
            test = load(test_file);
            control_tbl = unpack_table(control.res.table);
            label = cell(height(control_tbl),1);
            label(:) = {'control'};
            control_tbl.label = label;
            test_tbl = unpack_table(test.res.table);
            label = cell(height(test_tbl),1);
            label(:) = {'test'};
            test_tbl.label = label;
            obj.segment_tbl = vertcat(control_tbl,test_tbl);
            obj.n_px = control.res.n_px;
            obj.from_px = control.res.from_px;
            obj = obj.classify_opening;
        end
        function [new_obj_sub, new_obj_exc]  = subarea(obj,area_name)
            % Create new object with only sub area of the brain specified
            % as a string. example: new_obj = obj.subarea('hypothalamus');
            % Inputs:
            %   area_name (str)- string with the name of the ROI to be
            %       extracted
            area_idx = cellfun(@(x) contains(lower(x),lower(area_name)),...
                obj.segment_tbl.image_name);
            new_obj_sub = obj;
            new_obj_sub.segment_tbl(~area_idx,:) = [];
            new_obj_exc = obj;
            new_obj_exc.segment_tbl(area_idx,:) = [];
        end
        function writecsv(obj,ths,control_csv_filename,test_csv_filename,varargin)
            % save to csv files for statistical analysis in external SW
            % such as graphpad. the data will be saved into 2 seperate CSV
            % files, one for control and the other for MB+FUS. Each row of
            % the CSV represents a diameter group (e.g. 2-3) and each
            % column represents the intensity of the BBBD marker dye (EB)
            % around a single vessel segment belonging to the diameter
            % group.
            % 
            % Inputs:
            %   ths = diameter grouop edges
            %   control_csv_filename = path to csv file for control group
            %   test_csv_filename = path to csv file for test group
            % Name-Value pair arguements:
            %   GroupByFrame = boolean flag (default = 0) 
            P = inputParser();
            P.addOptional('GroupByFrame',0,@(x) ismember(x,[0,1]));
            P.parse(varargin{:});
            GroupByFrame = P.Results.GroupByFrame;
            if nargin < 1
                ths = 2:10;
            end
            if nargin < 3
                [control_csv_filename, control_csv_folder] = ...
                    uiputfile({'*.csv';'*.xlsx'},'Specify control csv file name');
                control_csv_filename = ...
                    fullfile(control_csv_folder,control_csv_filename);
                [test_csv_filename, test_csv_folder] = ...
                    uiputfile({'*.csv';'*.xlsx'},'Specify test csv file name');
                test_csv_filename = ...
                    fullfile(test_csv_folder, test_csv_filename);
            end
            control_idx = cellfun(@(x) strcmp(x,'control'),obj.segment_tbl.label);
            n_bins = numel(ths);
            control_discrete_cell = intogroups(...
                obj.segment_tbl(control_idx,:),ths,GroupByFrame);
            test_discrete_cell = intogroups(...
                obj.segment_tbl(~control_idx,:),ths,GroupByFrame);
            str_cell = cell(n_bins-1,1);
            ths = [0,ths];
            for i = 1:n_bins
                str_cell{i} = sprintf('%d - %d',ths(i),ths(i+1));
            end
            control_tbl = cell2table(...
                [str_cell,control_discrete_cell],...
                'VariableNames',{'Diameter','BBBD marker intensity'});
            test_tbl = cell2table(...
                [str_cell,test_discrete_cell],...
                'VariableNames',{'Diameter','BBBD marker intensity'});
            writetable(control_tbl, control_csv_filename);
            writetable(test_tbl, test_csv_filename);
        end
        function new_obj = keep_diameters(obj,lowLim,highLim)
           % Remove all vessel segments outside the given diameter range
           new_obj = obj;
           new_obj.segment_tbl(...
               (obj.segment_tbl.median_segment_diam_um < lowLim) |...
               (obj.segment_tbl.median_segment_diam_um > highLim),:) = [];
        end
        function new_obj = classify_opening(obj,ths,numstd)
            % Classify if a vessel was opened or not based on the intensity
            % of BBBD fluorescent marker in the perivascular area. the
            % threshold for opening is calculated as
            % mean{I_control}+numstd*SD{I_control}
            % The classification criteria is different for each diameter
            % group specified by ths.
            % Inputs:
            %    ths - diameter groups to be used for classification
            %           specified as a list of bin-edges, e.g. ths = [2,3] will
            %           create 3 diameter groups: [0-2],[2-3],[3-inf]
            %    numstd - number of standard deviations from the avarage
            %           BBBD marker intensity in the perivascular area IN
            %           THE CONTROL GROUP. This will be set as the
            %           threshold for BBBD opening.
            if nargin <3
                numstd = 2;
            end
            if nargin < 2
               ths = 2:10;
            end
            new_obj = obj;
            control_idx = cellfun(@(x) strcmp(x,'control'),...
                obj.segment_tbl.label);
            cont_ext_grouped = ...
                intogroups(obj.segment_tbl(control_idx,:),ths);
            cont_avs = cellfun(@(x) mean(x), cont_ext_grouped);
            cont_stds = cellfun(@(x) std(x), cont_ext_grouped);
            treat_th = cont_avs + numstd.*cont_stds; % Calculate the threshold
            % placeholder setup
            new_obj.segment_tbl.opening = ...
                zeros(height(new_obj.segment_tbl),1);
            ths = [0,ths];
            for i = 1:length(ths)-1
                new_obj.segment_tbl.opening(...
                    new_obj.segment_tbl.median_segment_diam_um >= ths(i) &...
                    new_obj.segment_tbl.median_segment_diam_um < ths(i+1) &...
                    ~control_idx &...
                    new_obj.segment_tbl.bbbd_marker_median >= treat_th(i)) = 1;
            end
        end
        function new_obj = remove_penetrating(obj, T_diam, T_len)
            % Function to remove vessels that have both small diameter and
            % short length (perpendicular to the image plane)
            % Inputs:
            %   T_diam = diameter threshold [um]
            %   T_len = length threshold [um]
            new_obj = obj;
            new_obj.segment_tbl(...
                new_obj.segment_tbl.median_segment_diam_um<T_diam &...
                new_obj.segment_tbl.len<T_len,:) = [];
            new_obj = new_obj.classify_opening;
        end
        %% Plotting functions
        function violinplot(obj,ths,varargin)
            % A violin plot of the BBBD marker intensity in the
            % perivascular area where the x-axis is the diameter group, the
            % y-axis is the intensity distribution. each violin is split 
            % to left (red) which represents control data, and right (blue)
            % which represents MB+FUS data.
            % Inputs:
            %   ths (optional)- array of diameters to be used as x-axis.
            %       vessels with diameter larger than ths(end) will not be
            %       presented. if "ths" is not specified a single violin 
            %       will be plotted for all vessels
            %   varargin - input arguements of the varargin function as
            %   presented here: github.com/bastibe/Violinplot-Matlab
            control_idx = cellfun(@(x) strcmp(x,'control'),...
                obj.segment_tbl.label);
            figure;
            if nargin < 2   % No diameter thresholds specified
                control_groups = ...
                    obj.segment_tbl.bbbd_marker_median(control_idx);
                test_groups = ...
                    obj.segment_tbl.bbbd_marker_median(~control_idx);
                control_groups = rmoutliers(control_groups);
                test_groups = rmoutliers(test_groups);
                Violin({control_groups},1,...
                    'HalfViolin','left','ViolinColor',{[1,0,0]},...
                    varargin{:});
                hold on;
                Violin({test_groups},1,...
                    'HalfViolin','right','ViolinColor',{[0,0,1]},...
                    varargin{:});
                xticks([]);
                hold off;
            else
                control_groups = ...
                    intogroups(obj.segment_tbl(control_idx,:),ths);
                test_groups = ...
                    intogroups(obj.segment_tbl(~control_idx,:),ths);
                % remove outliers
                control_groups = cellfun(@(x) rmoutliers(x),...
                    control_groups,'UniformOutput',false);
                test_groups = cellfun(@(x) rmoutliers(x),...
                    test_groups,'UniformOutput',false);
                for i = 1:numel(control_groups)
                    if ~isempty(control_groups{i})
                        Violin(control_groups(i),i,...
                            'HalfViolin','left','ViolinColor',{[1,0,0]});
                        hold on;
                    end
                end
                for i = 1:numel(test_groups)
                    if ~isempty(test_groups{i})
                        Violin(test_groups(i),i,...
                            'HalfViolin','right','ViolinColor',{[0,0,1]});
                        hold on;
                    end
                end
                xticks([1:numel(control_groups)]);
                xticklabels(generate_xticks(ths));
                xlabel('Diameter [um]');
                hold off;
            end
            ylabel('Median bbbd marker intensity in perivscular area [A.U.]');
            ax = gca;
            ch = get(ax,'Children');
            red_envalope = ch(end-1);
            blue_envalope = ch(7);
            legend([red_envalope,blue_envalope],{'control','MB + FUS'})
        end
        function barplot(obj,ths)
            % Bar plot of the intensity of BBBD marker in the perivascular
            % area, presented with significance stars.
            % Inputs:
            %   ths (optional)- array of diameters to be used as x-axis.
            %       vessels with diameter larger than ths(end) will not be
            %       presented.
            control_idx = cellfun(@(x) strcmp(x,'control'),...
                obj.segment_tbl.label);
            if nargin < 2   % No diameter thresholds specified
                control_groups = ...
                    obj.segment_tbl.bbbd_marker_median(control_idx);
                test_groups = ...
                    obj.segment_tbl.bbbd_marker_median(~control_idx);
                control_mu_median = [mean(control_groups);
                    std(control_groups)];
                test_mu_median = [mean(test_groups);
                    std(test_groups)];
                control_groups = rmoutliers(control_groups);
                test_groups = rmoutliers(test_groups);
                ths = 0;
            else
                control_groups = ...
                    intogroups(obj.segment_tbl(control_idx,:),ths);
                test_groups = ...
                    intogroups(obj.segment_tbl(~control_idx,:),ths);
                control_mu_median = cellfun(@(x) [mean(x);std(x)],...
                    control_groups,'UniformOutput',false);
                control_mu_median = [control_mu_median{:}];
                test_mu_median = cellfun(@(x) [mean(x);std(x)],...
                    test_groups,'UniformOutput',false);
                test_mu_median = [test_mu_median{:}];
                % remove outliers
                control_groups = cellfun(@(x) rmoutliers(x),...
                    control_groups,'UniformOutput',false);
                test_groups = cellfun(@(x) rmoutliers(x),...
                    test_groups,'UniformOutput',false);
            end
            if length(ths) == 1
                b1 = bar(0.75,control_mu_median(1,:),0.25,...
                    'FaceColor','#8c1515');
                hold on;
                b2 =bar(1.25,test_mu_median(1,:),0.25,...
                    'FaceColor','#09425A');
                errorbar(0.75,control_mu_median(1,:),...
                    control_mu_median(2,:),control_mu_median(2,:),...
                    'k', 'LineStyle','none');
                errorbar(1.25,test_mu_median(1,:),...
                    test_mu_median(2,:),test_mu_median(2,:),'k',...
                    'LineStyle','none'); 
                xticks([]);
                % Add significance stars of control vs test of same diameter
               [~,p] = ttest2(control_groups,test_groups);
               maxy = max([sum(control_mu_median),...
                   sum(test_mu_median)]);
               line([0.5,1.5],(maxy*1.05)*[1,1]);
               text(0.5,maxy*1.08,sigstars(p));
            else
                b1 = bar(0.75:2:(length(ths)*2),...
                    control_mu_median(1,:),0.25,...
                    'FaceColor','#8c1515');
                hold on;
                b2 =bar(1.25:2:(length(ths)*2),...
                    test_mu_median(1,:),0.25,...
                    'FaceColor','#09425A');
                errorbar(0.75:2:(length(ths)*2),control_mu_median(1,:),...
                    control_mu_median(2,:),control_mu_median(2,:),'k',...
                    'LineStyle','none');
                errorbar(1.25:2:(length(ths)*2),test_mu_median(1,:),...
                    test_mu_median(2,:),test_mu_median(2,:),'k',...
                    'LineStyle','none'); 
                xticks(1:2:(length(ths)*2));
                xticklabels(generate_xticks(ths));
                xlabel('Vessel diameter [um]');
                % Add significance stars of control vs test of same diameter
                for i = 1:length(ths)
                   [~,p] = ttest2(control_groups{i},test_groups{i});
                   maxy = max([sum(control_mu_median(:,i)),sum(test_mu_median(:,i))]);
                   x_cord = i-1;
                   line([0.5,1.5]+x_cord*2,(maxy*1.05)*[1,1]);
                   text(x_cord*2+0.5,maxy*1.08,sigstars(p));
                end
            end
            ylim([0,maxy*1.1]);
            legend([b1,b2],'control','MB + FUS');
            title(sprintf('[%d,%d] px',...
                obj.from_px,obj.from_px+obj.n_px));
            ylabel('Median BBBD marker intensity in perivascular area [A.U.]')     
        end
        function markerDistrebution(obj,ths,numstd,varargin)
            % Plot the distribution of bbbd marker intensity in the 
            % perivascular area. This function outputs 2 overlaying
            % histograms, red = control, blue = MB + FUS
            % Inputs:
            %   ths (optional)- array of diameters to be used as diameter ..
            %       groups.
            %   numstd (optional) - set and plot the threshold for BBBD 
            %       classification as described in "classify_opening"
            %       function
            % Name-Value pair arguements:
            %   'Mode' (optional) - display mode, can either be 'histogram'
            %       (default) or 'pdf' which displys the kernel density
            P = inputParser();
            P.addParameter('Mode','histogram',...
                @(x) sum(strcmp(x,{'histogram','pdf'}))==1);
            P.parse(varargin{:});
            showmode = strcmp(P.Results.Mode,'histogram');
            if nargin < 1
                ths = 2:10;
            end
            if nargin < 2
                numstd = 2;
            end
            control_idx = cellfun(@(x) strcmp(x,'control'),...
                obj.segment_tbl.label);
            figure;
            thresh_specified = (nargin>1);
            if thresh_specified
                control_groups = ...
                    intogroups(obj.segment_tbl(control_idx,:),ths);
                test_groups = ...
                    intogroups(obj.segment_tbl(~control_idx,:),ths);   
                len = numel(ths);
                ths = [0,ths];
            else
                control_groups =...
                    obj.segment_tbl.bbbd_marker_median(control_idx);
                test_groups = ...
                    obj.segment_tbl.bbbd_marker_median(~control_idx);
                len = 1;
            end
            for i = 1:len
                if thresh_specified
                    subplot(ceil(sqrt(len)),ceil(sqrt(len)),i);
                    cur_controls = control_groups{i};
                    cur_tests = test_groups{i};
                else
                    cur_controls = control_groups;
                    cur_tests = test_groups;
                end
                if ~isempty(cur_controls)
                    if showmode
                        a1 = histogram(cur_controls,100,'FaceColor','#8c1515',...
                            'Normalization','probability');
                        hold on;
                        a2 = histogram(cur_tests,100,'FaceColor','#09425A',...
                            'Normalization','probability');
                    else
                        [control_density, control_vals] =...
                            ksdensity(cur_controls);
                        [test_density, test_vals] =...
                            ksdensity(cur_tests);
                        a1 = area(control_vals,...
                            100*control_density./sum(control_density),...
                            'FaceColor','#8c1515');
                        a1.FaceAlpha = 0.5;
                        hold on;
                        a2 = area(test_vals,...
                            100*test_density./sum(test_density),...
                            'FaceColor','#09425A');
                        a2.FaceAlpha = 0.5;
                        [max_control,max_control_idx] =...
                            max(100*control_density./sum(control_density));
                        [max_test,max_test_idx] =...
                            max(100*test_density./sum(test_density));
                        plot([control_vals(max_control_idx),test_vals(max_test_idx)],...
                            (0.1 + max([max_control,max_test]))*ones(1,2),...
                            'Color',[0 0 0],'LineWidth',1);
                        all_measurements = vertcat(cur_controls,cur_tests);
                        all_label = ones(size(all_measurements));
                        all_label(1:length(cur_controls)) = 0;
                        p = anovan(all_measurements,all_label,'display','off');
                        text(...
                            mean([control_vals(max_control_idx),test_vals(max_test_idx)]),...
                            0.15 + max([max_control,max_test]),...
                            sigstars(p),...
                           'FontSize',12,...
                           'HorizontalAlignment','center');
                    end
                    xlabel('BBBD marker intensity [A.U.]');
                    ylabel('# Vessels [%]');
                    legend([a1,a2],{'control','MB + FUS'});
                    if isnumeric(numstd) && ~isempty(numstd)
                       l1 = xline(mean(cur_controls)+numstd*std(cur_controls),...
                           'LineWidth',2,'LineStyle','--');
                       legend([a1,a2,l1],{'control','MB + FUS',...
                           ['Control mean + ',num2str(numstd),' SDs']});
                    end
                end
                xlim([0,1]);
                ylim([0,5]);
                if thresh_specified
                    title(sprintf('%d-%d um diameter',ths(i),ths(i+1)));
                end
                hold off; 
            end
        end
        function diamHist(obj)
            % histogram of blood vessel diameter. This function outputs two
            % overlaying histograms, red = control vessels, blue = MB + FUS
            control_idx = cellfun(@(x) strcmp(x,'control'),...
                obj.segment_tbl.label);
            control_diams = obj.segment_tbl.median_segment_diam_um(control_idx);
            test_diams = obj.segment_tbl.median_segment_diam_um(~control_idx);
            histogram(control_diams,0.5:1:25.5,...
                'Normalization','probability','FaceColor','#007C92'); 
            hold on;
            histogram(test_diams,0.5:1:25.5,'Normalization','probability',...
                'FaceColor','#E98300');
            legend(['Control-',num2str(numel(control_diams)),' total segments']...
                ,['Treatment-',num2str(numel(test_diams)),' total segments']);
            xlabel('Vessel diameter [um]'); ylabel('% of total vessels');
            xticks(1:25)
            title('Blood vessel diameter histogram'); 
            p = anova1(obj.segment_tbl.median_segment_diam_um,~control_idx)
            xticklabels({'control','MB + FUS'});
            ylabel('Diameter [um]');
        end
        function perc_tbl = openedHist(obj,ths,varargin)
            % plot the histogram of the fraction of vesseles of each 
            % diameter that were classified as opened as per
            % Inputs:
            %   ths - array of diameters to be used as bin-edges for 
            %       grouping by diameter. default = [2:10]
            % Name-Value pair arguements:
            %   Errorbars - 'on' or 'off' (default)
            % Output:
            %   perc_tbl = opening percentage table by dimeter and frame            
            P = inputParser();
            P.addOptional('Errorbars','off',...
                @(x) sum(strcmp(x,{'on','off'})) == 1);
            P.parse(varargin{:})
            if nargin < 1
                ths = 2:10;
            end          
            control_idx = cellfun(@(x) strcmp(x,'control'),...
                obj.segment_tbl.label);
            test_frames = unique(obj.segment_tbl.image_name(~control_idx));
            ths = [0, ths];
            len_diam_groups = length(ths)-1;
            n_animals = numel(test_frames);
            perc = zeros(len_diam_groups,n_animals);
            vessel_count_per_brain = perc;
            for i = 1:len_diam_groups
                for j = 1:n_animals
                    in_group = ...
                        obj.segment_tbl.median_segment_diam_um >= ths(i) &...
                        obj.segment_tbl.median_segment_diam_um < ths(i+1) &...
                        strcmp(obj.segment_tbl.image_name,test_frames(j));
                    vessel_count_per_brain(i,j) = sum(in_group);
                    open_temp = in_group & obj.segment_tbl.opening;
                    perc(i,j) = 100*(sum(open_temp)/sum(in_group));
                end
            end
            bar(mean(perc,2,'omitnan'),0.5,'FaceColor','#009779');
            if strcmp(P.Results.Errorbars,'on')
                hold on;
                errorbar(1:len_diam_groups,mean(perc,2,'omitnan'),...
                            [],std(perc,0,2,'omitnan'),'k',...
                            'LineStyle','none');
            end
            xlabel('Blood vessel diameter [um]'); 
            xticklabels(generate_xticks(ths(2:end)));
            ylabel('Fraction of BBB opened vessels [%]');
            title('Opened vessel fraction as function of diameter'); 
            % return table
            perc_tbl = array2table(perc,'VariableNames',test_frames,...
                'RowNames',generate_xticks(ths(2:end)));
        end
    end
end

%% Helper functions
function detailed_tbl = unpack_table(tbl)
% convert a packed table to a new table where every row represent a single
% vessel segment
lens = cellfun(@(x) length(x),tbl.median_segment_diam_um);
cumlens = cumsum(lens);
imname = cell(sum(lens),1);
for i = 1:height(tbl)
    if i == 1
        start_idx = 1;
    else
        start_idx = cumlens(i-1)+1;
    end
    end_idx = cumlens(i);
    imname(start_idx:end_idx) = tbl.image_name(i);
end
detailed_tbl = table(imname,vertcat(tbl.median_segment_diam_um{:,:}),...
    vertcat(tbl.bbbd_marker_median{:,:}),vertcat(tbl.segment_len_um{:,:}),...
    'VariableNames',...
    {'image_name','median_segment_diam_um','bbbd_marker_median','len'}); 
end
function ticklabels = generate_xticks(ths)
% Generate a cell array of xticks based on ths. where each tick is a string
% of 'ths(i-1)-ths(i)' and for the 1st threshold it is '0-ths(1)'
tmp = cellfun(@(x) num2str(x),num2cell(ths),'UniformOutput',false);
tmp = [{'0'},tmp];
ticklabels = tmp(1:end-1);
for i = 2:numel(tmp)
    ticklabels(i-1) = ...
        cellstr([tmp{i-1},'-',tmp{i}]);
end
end
function str = sigstars(p)
% Function to calculate significance stars based on p val
if p<=10^-4
   str = '****';
elseif p <=10^-3
   str = '***';
elseif p <=10^-2
   str = '**';
elseif p <=0.05
   str = '*';
else
   str = 'ns';
end
end