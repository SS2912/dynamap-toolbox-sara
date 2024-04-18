function bigtable = MMTable(dir_data, fold, dir_info, name_info_table, varargin)

%% DATA ORGANIZATION IN TABLE FOR MIXED MODEL ANALYSIS
% Inputs: 
    % #1,2) dir_data, fold: folder with h2.mat results for each subject and
    % for a specific band (fold)
    % #3,4) dir_info, name_info_table: table with additional information for
    % each electrode contaning all subjects (ex: subj responder or not, anatomical region of contact, EI, ...)
    % OPTIONALS: 
    % a) 'sheets', 'subj_param'= sheet where to take info. If one sheet has info for
    % each electrode and another sheet has info for each subject, write it
    % like ["chan_sheet","subj_sheet"]
    % b) 'roi': region of interest inside which to calculate the strength/degrees values
    % options: "all", "ez", "pz", "ni", "ez+pz", "pz+ni". Default: all
    % c) 'thr_h2', 'thr_lag' : thr_h2=0 for strength, >0 for degrees. 
    % defaults: thr_h2=0, thr_lag=0
    % d) 'norm': do you wanna normalise the strength by the number of electrodes? yes=1/T, no=0/F
    
% Goal: create a table containing 1 row per each contact, specifying for
% each row the: subject, anat region, resp or not, type of epilepsy, etcetc
    
%% UPLOAD VALUES FILES
% Structure of data storage: 1 folder for each type of analysis,
% containing all matlab results from analysis

% organization of folder depends on what you wanna do. In this script, I
% want to calculate 1 value of strength (mean) for each channel for each
% subject and put it in the big table (bigT) so I put all subj's files
% together and a for cycle will let me calculate strengt and put it in the
% corresponding slot in bigT

% MAKE SURE YOU HAVE, file 1 = postTC2 (filepost(2).mat), file 2 = postTC1,  file 3 = preTC

% default values
thr_h2 = 0;  % thr_h2=0 for strength, >0 for degrees. Default: strength
thr_lag = 0; % value in ms, optional input. Default: 0 % TO CHECK???
sheets = "false";
bigT=[]; subj_struc= [];

% Optionals
    for ii = 1:2:nargin-4
        if strcmp('sheets', varargin{ii})
            clearvars sheets
            sheets = varargin{ii+1};    % options: none, one sheet, more sheets with different info as string vector
        elseif strcmp('roi', varargin{ii}) 
            roi_sel = varargin{ii+1};
        elseif strcmp('thr_h2', varargin{ii}) 
            thr_degreesH2 = varargin{ii+1};  
        elseif strcmp('thr_lag', varargin{ii}) 
            thr_lagH2 = varargin{ii+1};   
        elseif strcmp('subj_param', varargin{ii}) 
            subj_param = varargin{ii+1};  
        elseif strcmp('norm', varargin{ii}) %do you wanna normalise the strength by the number of electrodes? yes=1/T, no=0/F
            norm = varargin{ii+1};  
        end
    end
    

% 1. Upload additional information as excel sheet/table
whereistable = strcat(dir_info, '\', name_info_table);
if strcmp(sheets, "false")
    chan_info = readtable(whereistable); 
elseif length(sheets)==1
    chan_info = readtable(whereistable, 'Sheet', sheets);
elseif length(sheets)==2
    chan_info = readtable(whereistable, 'Sheet', sheets(1));
    subj_info = readtable(whereistable, 'Sheet', sheets(2));
end

% 2. set directory for h2 files to analyse
folder= strcat(dir_data,'\',fold);
cd(folder)
myfiles = dir('*.mat');
% files are organized as: 1= postTC1; 2= postTC2; 3=preTC


% 3. Calculate the strength/degrees for each channel and each patient and
% each channel

% Option --> Calculate between selected channels (in the "roi" input option): you can
% calculate stre/deg only between EZ or PZ or NI, between all of them
% or pooling together two of them (e.g. EZ+PZ)

switch roi_sel
    case "ez"
        good_labels= "EZ";
    case "pz"
        good_labels= "PZ";
    case "ni"
        good_labels= "NI";
    case "ez+pz"
        good_labels= '(EZ|PZ)';
    case "pz+ni"
        good_labels= '(NI|PZ)';
    otherwise
        good_labels = "all";
end
    

for index=1:3:length(myfiles)
    datapostB = load(myfiles(index).name); % in old results (feb-march 2021), postA and postB are reversed in the order. 
    datapostA = load(myfiles(index+1).name);
    datapre = load(myfiles(index+2).name);
    subj = string(extractBetween(myfiles(index).name, 'sub-', '_ses'));
    channels = string(datapre.electrode_names);
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
    % 3.2. Reduce matrix of h2 and lag in case "roi" is different from "all"
    if ~strcmp("all", good_labels)
        H={datapre, datapostA, datapostB};
        for dd= 1:3
            idx = ~cellfun(@isempty,regexp(subj_struc(:,3),good_labels));
            H{dd}.aw_h2 = H{dd}.aw_h2(idx,idx,:);
            H{dd}.aw_lag = H{dd}.aw_lag(idx,idx,:);
            H{dd}.electrode_names = H{dd}.electrode_names(idx);
        end
        datapre = H{1}; datapostA = H{2}; datapostB = H{3};
        channels = string(datapre.electrode_names); %update the selected channels
        subj_struc = subj_struc(idx,:);
    end
    
    % add additional information (from subj_info and the parameters specified as input) to the subject matrix
    other_info = subj_info(subj_info.BIDSID == subj, subj_param);
    other_info = string(table2cell(other_info));
    subj_struc = [subj_struc, repmat(other_info, length(channels),1)];
    
    % 3.3 calculate the in, out, tot strengths/degrees using the func countlinks in graphcompare
    if thr_h2==0 && norm==1
        
        [linkspre,~]=ins_countlinks(datapre,thr_h2);        % linkspre contains 1 row per channel, 1 column per window of h2
        OUTpre = median(linkspre.outstrength_norm,2);
        TOTpre = median(linkspre.totstrength_norm,2);
        %post1
        [linkspostA,~]=ins_countlinks(datapostA,thr_h2);
        OUTpostA = median(linkspostA.outstrength_norm,2);
        TOTpostA = median(linkspostA.totstrength_norm,2);
        %post2
        [linkspostB,~]=ins_countlinks(datapostB,thr_h2);
        OUTpostB = median(linkspostB.outstrength_norm,2);
        TOTpostB = median(linkspostB.totstrength_norm,2);
        
        for chan = 1: size(linkspre.totstrength,1)  
            clearvars statsa statsb
            [~,~,statsa] = ranksum(linkspostA.totstrength_norm(chan,:), linkspre.totstrength_norm(chan,:)); %zvalue of norm is same of not norm (I also checked that)
            ztota(chan) = statsa.zval;
            [~,~,statsb] = ranksum( linkspostB.totstrength_norm(chan,:), linkspre.totstrength_norm(chan,:));
            ztotb(chan) = statsb.zval;

        end
        
        TOTzA = ztota';
        TOTzB = ztotb';
       clearvars ztota ztotb
        
    else if thr_h2==0 && norm==0
            
        [linkspre,~]=ins_countlinks(datapre,thr_h2);
        OUTpre = median(linkspre.outstrength,2);
        TOTpre = median(linkspre.totstrength,2);
        %post1
        [linkspostA,~]=ins_countlinks(datapostA,thr_h2);
        OUTpostA = median(linkspostA.outstrength,2);
        TOTpostA = median(linkspostA.totstrength,2);
        %post2
        [linkspostB,~]=ins_countlinks(datapostB,thr_h2);
        OUTpostB = median(linkspostB.outstrength,2);
        TOTpostB = median(linkspostB.totstrength,2);
       
        for chan = 1: size(linkspre.totstrength,1)  
            clearvars statsa statsb
            [~,~,statsa] = ranksum(linkspostA.totstrength(chan,:), linkspre.totstrength(chan,:));
            ztota(chan) = statsa.zval;
            [~,~,statsb] = ranksum(linkspostB.totstrength(chan,:), linkspre.totstrength(chan,:));
            ztotb(chan) = statsb.zval;

        end
        
        TOTzA = ztota';
        TOTzB = ztotb';
       clearvars ztota ztotb
        
    else 
        [linkspre,~]=ins_countlinks(datapre,thr_h2);
        OUTpre = median(linkspre.outdegree,2);
        TOTpre = median(linkspre.totdegree,2);
        %post1
        [linkspostA,~]=ins_countlinks(datapostA,thr_h2);
        OUTpostA = median(linkspostA.outdegree,2);
        TOTpostA = median(linkspostA.totdegree,2);
        %post2
        [linkspostB,~]=ins_countlinks(datapostB,thr_h2);
        OUTpostB = median(linkspostB.outdegree,2);
        TOTpostB = median(linkspostB.totdegree,2);
        
        for chan = 1: size(linkspre.totstrength,1)  
            clearvars statsa statsb
            [~,~,statsa] = ranksum(linkspre.totdegree(chan,:), linkspostA.totdegree(chan,:));
            ztota(chan) = statsa.zval;
            [~,~,statsb] = ranksum(linkspre.totdegree(chan,:), linkspostB.totdegree(chan,:));
            ztotb(chan) = statsb.zval;

        end
        
        TOTzA = ztota';
        TOTzB = ztotb';
       clearvars ztota ztotb

    end
    end
    % concatenate the subj info table with the h2 values
    hvalues = string([OUTpre; TOTpre; OUTpostA; TOTpostA; OUTpostB; TOTpostB; TOTzA; TOTzB]);
    hvalues(:,2) =  repelem(["OUTpre", "TOTpre", "OUTpostA", "TOTpostA",...
        "OUTpostB", "TOTpostB", "TOTzvalA", "TOTzvalB"]', length(channels),1);
    
    temp = [repmat(subj_struc,8,1) hvalues];
    bigT = [bigT; temp];

    
    clearvars hvalues datapre datapostA datapostB subj channels check linkspre linkspostA linkspostB subchan OUTpre TOTpre OUTpostA TOTpostA OUTpostB TOTpostB H subj_struc TOTzvalA TOTzvalB
end

bigtable = array2table(bigT,...
    'VariableNames',{'subj','channel','roi','thermo','responder', ...
    subj_param{:}, 'h2_value', 'h2_type'});

% name = strcat(dir_data,'\',fold,'_',roi_sel,'.csv');
% writetable(bigtable,name)  

end

