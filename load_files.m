function [allfiles]=load_files(string, samestructure)

%% Loads only files with chosen string in the file name

myfile = dir(string);  %reads the properties of all files containing the chosen string in their name (and with .mat extension if specified)

%IN case there are multiple files to be loaded but they have the same
%variable names, you have to put them in a structure so to avoid
%overwriting of previous variables --> choose samestructure as 1
if samestructure ==1
    n= length(myfile);
        for i=1:n         %sometimes I have diffenent time periods before and/or after stimulation
            allfiles(i)=load(myfile(i).name);   %loads the structures 'result' and 'detection' in the allfiles structure
        end
    
else       %if variables have different names, you can load them in separate variables
    n= length(myfile);
        for i=1:n         %sometimes I have diffenent time periods before and/or after stimulation
            load(myfile(i).name);   %loads the structures 'result' and 'detection' in the ALL_before structure
        end
end
