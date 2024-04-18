function zplot(data, fmax, window, sortby, movavg, colorlim)
% computes and plot the zscore values for the difference/ ratio matrices of
% spectral density for each channel. Plots a figure with 1 subplot for each
% patient with code and responder dets. If movavg is 1, plots the zscore
% values computed with moving mean along each channel

% author: Sara S. thermocoagulation project, summer 2021
if ~exist('colorlim')
    colorlim = [-4 4];
end

for i = 1:size(data,2)
    
    if ~movavg
        data(i).zscorediffA = zscore(data(i).psddiffA);
        data(i).zscorediffB = zscore(data(i).psddiffB);
    
    else
        data(i).zscorediffA = movmean(zscore(data(i).psddiffA), 5);  %no need to put dimension = 2 cause cols=channels -- compute by row
        data(i).zscorediffB = movmean(zscore(data(i).psddiffB), 5);
    end
    
    a = ceil(sqrt(size(data,2)));
    b = floor(sqrt(size(data,2)));
    if  a*b < size(data,2)
        b = a;
    end        
    
    switch sortby
        case 'TC'     %sort by TC channel
        [~, indexes] = sort(data(i).chaninfo(:,4), 'descend');
        case 'ROI'    %or sort by involved/non involved
        [~, indexes] = sort(extractAfter(data(i).chaninfo(:,3),1), 'descend');
    end
    
    data(i).zscorediffB = data(i).zscorediffB(:,indexes);
    
    subplot(a, b, i);
    imagesc([data(i).zscorediffB]')
    xticks(0:fmax*window/10:fmax*window);
    xticklabels(0:fmax/10:fmax);
    yticks(1:1:size(data(i).zscorediffB, 2));
    yticklabels(strcat(data(i).chaninfo(indexes,2), '-', data(i).chaninfo(indexes,3), '-', data(i).chaninfo(indexes,4)));
    
    colorbar
    caxis(colorlim)
    title(strcat(data(i).code, '-', data(i).resp));
    ylabel('channel label');
    xlabel('frequency (Hz)');
end
end
