function [detections]=select_markers(results,marker)

%% Analyze DELPHOS results with sliding window
% 
% Syntax:  
%    [detected_markers]=select_markers(results,marker)
%
% Inputs:
%   results    - matlab file obtained with delphos
%   marker     - 'Spike', 'Ripple, 'Fast Ripple', 'Cross Rate', 'HFO'
% 
% Outputs:
%   detected_markers      - All the detections (timepoints) for each
%                        channel and for the selected marker type
%
% See also:

% Original author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% New version made by Sara specific for marker selection

%
if strcmp(marker,'Cross Rate')
    detections.spike = select_markers(results,'Spike');
    detections.ripple = select_markers(results,'Ripple');
    detections.fastripple = select_markers(results,'Fast Ripple');
    return
end
   
if strcmp(marker,'HFO')
    detections_aux.ripple = select_markers(results,'Ripple');
    detections_aux.fastripple = select_markers(results,'Fast Ripple');
    detections=cell(length(results.labels),1);
    for i=1:length(detections_aux.ripple)
        detections{i} = sort([detections_aux.ripple{i} detections_aux.fastripple{i}]);
    end
    return
end

detections=cell(length(results.labels),1);  %initializes the cell of spikes with the length equal to the number of channels present


% same cycle I made for extracting all spikes from the structure
for iter=1:length(results.markers)   %each signle marker detected
    if strcmp(results.markers(iter).label,marker) == 1
        for iter2=1:length(results.labels)
            if strcmp(results.markers(iter).channels, results.labels(iter2)) == 1
                detections{iter2} = [detections{iter2} results.markers(iter).position];
            end
        end
    end
end

%%if only selection of struct and no reorganisation:
% marker = "Spike";
% 
% j=1; a=ALL_before(1).results;
% for i=1:length(a.markers)
%     clear marker_name;
%     marker_name =a.markers(i).label;
%     if marker_name == marker
%         selectionA_before(j)=a.markers(i); j=j+1;
%     end
% end

