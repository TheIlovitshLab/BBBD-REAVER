function eb_grouped = intogroups(tbl,edges,GroupByFrame)
% Function to devide vessles into groups based on diameter bin edges
% Input arguements:
%     tbl = table object containing a mean_eb_diam column and a avg_eb_ext
%       column
%     edges = vector of bin edges for diameters
%     GroupByFrame = boolean flag, wether to average wach frame by itself 
% Output arguemrnts:
%     eb_grouped = cell vector of length numel(th)+1, where each cell
%     contains the eb values of segments belogning to the specific bin
if nargin < 3
    GroupByFrame = 0;
end
eb_grouped = cell(length(edges),1);
edges = [0,edges];
if ~GroupByFrame
    for i = 1:length(edges)-1
        eb_grouped{i} = ...
            tbl.median_red(edges(i)<=tbl.median_segment_diam_um &...
            tbl.median_segment_diam_um<edges(i+1));
    end
else
    frame_names = unique(tbl.image_name);
    for i = 1:length(edges)-1
        in_group = tbl(edges(i)<=tbl.median_segment_diam_um &...
                tbl.median_segment_diam_um<edges(i+1),:);
        eb_by_frame = zeros(numel(frame_names),1);
        for f = 1:numel(frame_names)
            eb_by_frame(f) = ...
                mean(...
                in_group.median_red(...
                strcmp(in_group.image_name,frame_names(f))...
                )...
                ,'omitnan');            
        end
        eb_by_frame(isnan(eb_by_frame)) = [];
        eb_grouped{i} = eb_by_frame;
    end
end