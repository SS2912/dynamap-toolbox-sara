function [h2_mean, h2_median] = FCmatrix_no_rep(h2mat, h2lag, anat_subj, chan_subj, chan_h2)
%
% Calculates mean or median across each structure to obtain only 1 value per structure 
% Syntax:  
%    FCtable = FCmatrix_no_rep(h2mat, h2lag, anat)
%
% Inputs:
%   h2mat       - aw.h2 matrix of interest (3D)
%   h2lag       - aw.lag matric of interest (3D)
%   anat        - column vector of all regions of the current subject
%   chan_subj   - column vector of all channels of the current subject
%   chan_h2     - vector of all channels from h2 .mat results (h2.electrode_names)
% 
% Varargin:
%   
% Output:
%   h2_mean     - struct with new h2 (h2_mean.h2) and lag matrix (h2_mean.lag)
%   h2_median   - struct with new h2 (h2_median.h2) and lag matrix (h2_median.lag)
%
% Required functions: 
%   
% Authors: Sara Simula (original: March 2024. Last version: )


%% Matrix symmetry
% necessary for median
% optional for mean: I compared both symm and non-symm results and they are
% the same but the scale changes a bit, with always higher h2 values for
% symmetrized matrix, but consistent in all patients. 
% Visible in \\dynaserv\Shared\Maria\Joy_project\Sara_analyses\1.rawdata\delta
for matrix_iter = 1 : size(h2mat,3)
    Wmatrix = h2mat(:,:,matrix_iter);
    for i = 1 : length(Wmatrix)
        for j = 1 : length(Wmatrix)
            if abs(Wmatrix(j,i))>=abs(Wmatrix(i,j))
                Wmatrix(i,j) = Wmatrix(j,i);
            end
        end
    end
    h2mat_symm(:,:,matrix_iter) = abs(Wmatrix); 
end

%% Initialization of vars and median/mean
A_mean = h2mat_symm;
A_median = h2mat_symm;
A_lag_mean = h2lag; %asymmetric for now but I dont actually care abt lag
A_lag_med = h2lag;
h2_mean = struct();
h2_median = struct();

chan_h2 = upper(chan_h2);
chan_subj = upper(chan_subj);
anat_h2 = strings([length(chan_h2),1]);
idx = zeros(length(chan_subj), 1);

for chan=1:length(chan_h2)
    idx = idx + contains(chan_subj, chan_h2(chan));
end

anat_h2 = anat_subj(logical(idx));
[uniqueStr, ~] = unique(anat_h2);
idx_norep = 0;
nbrwin = size(A_mean, 3);
new_anat =[];

for i = 1:length(uniqueStr)
    idx_norep = idx_norep +1;
    rep = contains(anat_subj, uniqueStr(i));
    rep_chan = find(contains(chan_h2, chan_subj(rep)));

    if length(rep_chan) > 1
        new_size = size(A_mean,1) - length(rep_chan) + 1;
        idx_keep = ~ismember(1:size(A_mean,1), rep_chan);
        new_chan = strings(1, new_size);
        new_chan(rep_chan(1)) = chan_h2(rep_chan(1));
        new_chan(~ismember(1:length(new_chan), rep_chan(1))) = chan_h2(idx_keep);
        new_anat = [new_anat, uniqueStr(i)];

        %% median
        B_median = zeros(new_size, size(A_median,1), nbrwin);
        B_median(idx_norep, :, :) = median(A_median(rep_chan, :, :)); %obtain row of median across same region with original columns
        B_median(~ismember(1:size(B_median,1), idx_norep), :, :) = A_median(idx_keep, :, :);

        B_lag_med = zeros(new_size, size(A_lag_med,1), nbrwin);
        B_lag_med(idx_norep, :, :) = median(A_lag_med(rep_chan, :, :)); %obtain row of median across same region with original columns
        B_lag_med(~ismember(1:size(B_lag_med,1), idx_norep), :, :) = A_lag_med(idx_keep, :, :);               

        %median over columns 
        C_median = zeros(new_size, new_size, nbrwin);
        C_median(:, idx_norep, :) = median(B_median(:, rep_chan, :), 2); %obtain column of median across same region with original columns
        C_median(:, ~ismember(1:size(C_median,1), idx_norep), :) = B_median(:, idx_keep, :);

        C_lag_med = zeros(new_size, new_size, nbrwin);
        C_lag_med(:, idx_norep, :) = median(B_lag_med(:, rep_chan, :), 2); %obtain column of median across same region with original columns
        C_lag_med(:, ~ismember(1:size(C_lag_med,1), idx_norep), :) = B_lag_med(:, idx_keep, :);

        %% mean
        %mean over rows 
        B_mean = zeros(new_size, size(A_mean,1), nbrwin);
        B_mean(idx_norep, :, :) = mean(A_mean(rep_chan, :, :)); %obtain row of mean across same region with original columns
        B_mean(~ismember(1:size(B_mean,1), idx_norep), :, :) = A_mean(idx_keep, :, :);

        B_lag_mean = zeros(new_size, size(A_lag_mean,1), nbrwin);
        B_lag_mean(idx_norep, :, :) = mean(A_lag_mean(rep_chan, :, :)); %obtain row of mean across same region with original columns
        B_lag_mean(~ismember(1:size(B_lag_mean,1), idx_norep), :, :) = A_lag_mean(idx_keep, :, :);               

        %mean over columns 
        C_mean = zeros(new_size, new_size, nbrwin);
        C_mean(:, idx_norep, :) = mean(B_mean(:, rep_chan, :), 2); %obtain column of mean across same region with original columns
        C_mean(:, ~ismember(1:size(C_mean,1), idx_norep), :) = B_mean(:, idx_keep, :);

        C_lag_mean = zeros(new_size, new_size, nbrwin);
        C_lag_mean(:, idx_norep, :) = mean(B_lag_mean(:, rep_chan, :), 2); %obtain column of mean across same region with original columns
        C_lag_mean(:, ~ismember(1:size(C_lag_mean,1), idx_norep), :) = B_lag_mean(:, idx_keep, :);
        
        
        % Set diagonal elements of each slice to 0
        for k = 1:nbrwin
            for diag = 1:size(C_median,1)
                C_median(diag, diag, k) = 0;
                C_lag_med(diag, diag, k) = 0;

                C_mean(diag, diag, k) = 0;
                C_lag_mean(diag, diag, k) = 0;
            end
        end   
        

        clearvars A_mean A_lag_mean A_median A_lag_med chan_h2
        A_mean = C_mean;
        A_lag_mean = C_lag_mean;
        A_median = C_median;
        A_lag_med = C_lag_med;
        chan_h2 = new_chan;
    
    end

    h2_median.h2 = A_median;
    h2_median.lag = A_lag_med;
    h2_median.chan = chan_h2;

    h2_mean.h2 = A_mean;
    h2_mean.lag = A_lag_mean;
    h2_mean.chan = chan_h2;
    
end

end