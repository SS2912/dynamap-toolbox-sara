function varargout = load_conn_matrix(dir_data, fold, matrix_names, path_info, varargin)

%%% DATA ORGANIZATION IN TABLE FOR MIXED MODEL ANALYSIS %%%
% Inputs: 
    % #1,2) dir_data, fold: folder with h2.mat results for each subject and
    % for a specific band (fold)
    % #3) matrix_names: which matrices you are uploading (ex: {"pre",
    % "post", "postb"} or {"pre", "during", "post"})
    % #4) 'path_info': table with additional information for
    % each electrode contaning all subjects (ex: subj responder or not, anatomical region of contact, EI, ...)
    % OPTIONALS: 
    % 'sheets', 'subj_param'= sheet where to take info. If one sheet has info for
    % each electrode and another sheet has info for each subject, write it
    % like ["chan_sheet","subj_sheet"]
    
% Goal: load the connectivity matrix for each patient, make h2 matrix symmetrical and 
% create a table containing 1 row per each contact, specifying for each row the: 
% subject, anat region, resp or not, type of epilepsy, etcetc
    
%% UPLOAD VALUES FILES
% Structure of data storage: 1 folder for each type of analysis,
% containing all matlab results from analysis. matrix_names tells the
% function how many data per patient to take


% default values
sheets = "false";
bigT=[]; subj_struc= [];
roi_sel = "all";

% Optionals
    for ii = 1:2:nargin-4 % 4 mandatory input vars (diardata, fold, matrix_names, path_info)
        if strcmp('sheets', varargin{ii})
            clearvars sheets
            sheets = varargin{ii+1};    % options: none, one sheet, more sheets with different info as string vector
        elseif strcmp('subj_param', varargin{ii}) 
            subj_param = varargin{ii+1};  
        elseif strcmp('path_info', varargin{ii})
            path_info = varargin{ii};
        end
    end

% 1. set directory for h2 files to analyse
folder= strcat(dir_data,'\',fold);
cd(folder)
myfiles = dir('*.mat');
% check that files are organised as in matrix_names
fprintf("Check that matrix_names correspond to this order: \n %s \n %s \n %s \n %s \n %s", string({myfiles(1:5).name}'))
answer = questdlg('Are matrix_names correct?', ...
	'Check organisation',...
    'Yes', 'No, stop running', 'Yes');
% Handle response
switch answer
    case 'Yes'
        disp([answer 'loading matrices...'])
    case 'No, stop running'
        return
end


% 2. Upload additional information as excel sheet/table
if strcmp(sheets, "false")
    chan_info = readtable(path_info); 
elseif length(sheets)==1
    chan_info = readtable(path_info, 'Sheet', sheets);
elseif length(sheets)==2
    chan_info = readtable(path_info, 'Sheet', sheets(1));
    subj_info = readtable(path_info, 'Sheet', sheets(2));
end

% 3. Create cell with aw_h2 matrices for each patient (row) and each condition (column)
idx = 0;   
conn_mtx = cell(length(myfiles)/length(matrix_names), length(matrix_names)+2); % +2 cause there is also subj + channel info in 2 additional columns

for sub_idx=1:length(matrix_names):length(myfiles)
    idx = idx+1;
    % for each condition(file, e.g. pre or post or during etcetc) load matrix and make it symmetric by taking max of the h2 connections
    for file = 1:length(matrix_names)  
        H = [];
        H = load(myfiles(sub_idx+file-1).name);
        temp = [];
        for z=1:size(H.aw_h2, 3)
            maxH = [];
            maxH = max(triu(H.aw_h2(:,:,z)), tril(H.aw_h2(:,:,z))');
            temp(:,:,z) = maxH + maxH'; %NB: to change with tru(maxH,+1)' if diagonal =/= 0, but DIAGONAL SHOUDL ALWAYS be 0!!
        end
        conn_mtx{idx,file} = temp;
    end 
    
    subj = string(extractBetween(myfiles(sub_idx).name, 'sub-', '_ses'));
    conn_mtx{idx,length(matrix_names)+1} = subj;
    channels = string(H.electrode_names);
    conn_mtx{idx, length(matrix_names)+2} = channels';
    conn_mtx{idx, length(matrix_names)+3} = matrix_names;
    subj_struc = [];
    
    % 3.1. Associate each channel to its ROI (&/or other info) from the info_table
    subchan = chan_info(chan_info.SubjBIDS == subj, :); %extract subtable of subject's data from the excel sheet
    % for each channel in 'channels', add info present on the excel table
    for i=1:length(channels)
       check=string(channels(i));
       if ~isempty(string(subchan.ROI(subchan.Channel == check)))
           subj_struc = [subj_struc; subj, check, string(subchan.ROI(subchan.Channel == check)), ...
               string(subchan.Thermo(subchan.Channel == check)), string(subchan.Responder(subchan.Channel == check))];
       else
           subj_struc = [subj_struc; subj, check, "NaN", "", ""];
       end    
    end
    
    
    % add additional information (from subj_info and the parameters specified as input) to the subject matrix
    other_info = subj_info(subj_info.BIDSID == subj, subj_param);
    other_info = string(table2cell(other_info));
    subj_struc = [subj_struc, repmat(other_info, length(channels),1)];
    bigT = [bigT; subj_struc];
    
end
    
    Varnames = ["subj", "channel", "roi", "thermo", "resp", subj_param];
    T = array2table(bigT, 'VariableNames', Varnames);
    
    varargout{1} = conn_mtx;
    varargout{2} = T;
   
end
    
