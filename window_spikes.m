
function [counts_window, mean_window, spikes2] = window_spikes(results, windowsize, varargin)

%% Analyze DELPHOS results with sliding window
% 
% Syntax:   
%    [mean_window,spikes]=window_spikes(results,windowsize)
%
% Inputs:
%   results    - matlab file obtained with DELPHOS
%   windowsize - length of the window in seconds
%   varargin: 
        % event   - type of event to count (e.g. 'Spike', 'Ripple', 'Gamma', etc. - default: 'Spike')
        % overlap - desired overlap for calculation (default = 0; should be 0 for counting spikes and hfos)
% 
% Outputs:
%   mean_window - [Nch x Nw] Mean rate at each window and for each channel
%   spikes2     - All the detections (timepoints) for each channel
%
% See also:

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Jan. 2020; Last revision by VLM: 22-Jan-2020
% Changes by SS, 12-Oct-2021 (added counts_window + changes 'spikes' into 'events' so that calculation of HFOs events is possible)

% default values 
event   = 'Spike';
overlap = 0;

for i= 1:nargin-2
    switch varargin{i}
        case 'event'
            event = varargin{i+1};
        case 'overlap'
            overlap = varargin{i+1};
    end
end

spikes=cell(length(results.labels),1);
spikes2=[];

for iter=1:length(results.markers)
    if strcmp(results.markers(iter).label, event) == 1
        for iter2=1:length(results.labels)
            if strcmp(results.markers(iter).channels, results.labels(iter2)) == 1
                spikes{iter2} = [spikes{iter2} results.markers(iter).position];
                spikes2 = [spikes2 results.markers(iter).position];
            end
        end
    end
end

step = round(windowsize*(1-overlap));
init_time=results.cfg.start;
Nw = round((results.cfg.duration-windowsize)/step+1);

for w=1:Nw
    init = init_time + (w-1)*step;
    endt = init + windowsize;
    
    for ch=1:length(results.labels)
        vector_of_spikesperwindow = spikes{ch}>init & spikes{ch}<endt;
        counts_window(ch,w) = sum(vector_of_spikesperwindow);
        mean_window(ch,w) = sum(vector_of_spikesperwindow)/windowsize;
    end
end

