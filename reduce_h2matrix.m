function [h2_selection, subj_info_sel] = reduce_h2matrix(h2_raw, roi, subj_info, roi_col)
%
% Reduce h2 matrix to contain only the regions selected in "roi"
% Syntax:  
%    [h2_selection, subj_info_sel] = reduce_h2matrix(h2_raw, roi, subj_info, roi_col)
%
% Inputs:
%   h2_raw      - struct containing aw_h2, aw_lag and electrode_names
%   roi         - string with desired selection (must be in a column of subj_info. examples: "EZ", "EZPZ", "NI", "involved")
%   subj_info   - table with all info of that subject (one row per channel)
%   roi_col     - index (number) of the subj_info's column where the label is stored
% 
% Varargin:
%   
% Output:
%   h2_selection    - updated h2 structure with selected channels only
%   subj_info_sel   - updated subj_info table with selected channels only
%
% Required functions: 
%   
% Authors: Sara Simula (original: March 2024. Last version: )

h2_selection = struct();

if strcmp(roi, "EZPZ")
    clearvars roi
    roi = '[EP.]Z';
end

idx = ~cellfun(@isempty, regexp(subj_info(:, roi_col), roi));
h2_selection.aw_h2 = h2_raw.aw_h2(idx,idx,:);
h2_selection.aw_lag = h2_raw.aw_lag(idx,idx,:);
h2_selection.electrode_names = h2_raw.electrode_names(idx);

subj_info_sel = subj_info(idx,:);

end
