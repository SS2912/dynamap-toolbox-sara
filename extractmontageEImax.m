function extractmontageEImax(numfiles, input, extension)

% extractmontageEImax creates a montage .mtg file, taking the 1 electrode with highest EI
% from each structure(from VEP atlas). This file can then be uploaded in AnyWave
% to reduce the montage for very heavy computations (e.g. h2 index)
% numfiles: 'single' for a single file or 'folder' for multiple files sored
% in one single folder
% input: name of the file (+extension) if 'single'; path of the folder if
% 'folder'
% extension: 'tsv' or 'xlsx' (more can be added if needed)


n=4; % specify the number of columns you wanna take 
        % n=4 takes only channels, EI, structure and zone (ROI)
        
        
if strcmp(numfiles, 'folder') 
    cd (input)
    switch extension
        case 'tsv'
            myfile = dir('*.tsv');  %reads the properties of all files containing the chosen string in their name (and with .mat extension if specified)
        case 'xlsx'
            myfile = dir('*.xlsx');
        otherwise disp ('wrong or missing extension: choose tsv or excel')
    end
    
else
    if strcmp(numfiles, 'single')
        myfile = {};
        myfile(1).name = input;
    end
end



switch extension 
    case 'tsv'     %if tsv file: this function needs the "tsvread" function     
     
    for i=1:length(myfile)        %sometimes I have diffenent time periods before and/or after stimulation
            [~, ~, data] = tsvread(myfile(i).name);
            VarNames= matlab.lang.makeValidName(string(data(1,1:n)));
            data= data(2:end, 1:n); 
            T= cell2table(data, 'VariableNames', VarNames);
            
            %% extract channels with maxEI, 1 for each region
            Tsort=sortrows(T, 'EIMax', 'descend');
            if strcmpi(Tsort(:,'Structure'), {'WM','White Matter', 'whitematter'}) ~= [0 0 0]
            [~, ia]= unique(Tsort(:,'Structure'), 'stable');
            MontageEImax= Tsort(ia, :); 

            channels= table2array(MontageEImax(:,'Channel'));%first line is column names so i take form row2
            
            %create file name for montage with subject code
            subject = extractBefore(myfile(i).name, ".tsv");
            namemtg = sprintf('mtgEImax_%s.mtg',subject);
            
            %% create a mtg file and write each selected bipolar elec in it
            fid = fopen(fullfile(cd,namemtg),'wt');
            % if fid
            fprintf(fid, '<!DOCTYPE AnyWaveMontage>\n<Montage>\n');
            for j = 1:size(channels,1)
                labels = strsplit(channels{j},'-');
                fprintf(fid, '\t<Channel name="%s">\n\t\t<type>SEEG</type>\n\t\t<reference>%s</reference>\n\t\t<color>black</color>\n\t</Channel>\n',...
                        labels{1},labels{2});
            %end
            end
            fprintf(fid, '</Montage>');
            fclose(fid);
        end
    case 'xlsx'
        %% same thing as tsv but simpler reading of xlsx file:
        for i=1:length(myfile)  
            T=readtable(myfile(i).name);
            
            Tsort=sortrows(T, 'EIMax', 'descend');
            [~, ia]= unique(Tsort(:,'Structure'), 'stable');
            MontageEImax= Tsort(ia, :); 
            
            channels = table2array(MontageEImax(:,'Channel')); %first line is column names so i take form row2
            
            %create file name for montage with subject code
            subject = extractBefore(myfile(i).name, ".tsv");
            namemtg = sprintf('mtgEImax_%s.mtg',subject);
            
            fid = fopen(fullfile(cd,namemtg),'wt');
            % if fid
            fprintf(fid, '<!DOCTYPE AnyWaveMontage>\n<Montage>\n');
            for j = 1:size(channels,1)
                labels = strsplit(channels{j},'-');
                fprintf(fid, '\t<Channel name="%s">\n\t\t<type>SEEG</type>\n\t\t<reference>%s</reference>\n\t\t<color>black</color>\n\t</Channel>\n',...
                        labels{1},labels{2});
            %end
            end
            fprintf(fid, '</Montage>');
            fclose(fid);
        end
end
end



