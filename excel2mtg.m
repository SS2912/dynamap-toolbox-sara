function excel2mtg (filename)

% trasnforms an excel file with channel names into an Anywave .mtg file

subject = extractBefore(filename, ".xlsx");
namemtg = sprintf('mtgEImax_%s.mtg',subject);       

[num,txt,raw] = xlsread(filename); 

channels= raw(2:end,1); %first line is column names

fid = fopen(fullfile(cd,namemtg),'wt');
% if fid
fprintf(fid, '<!DOCTYPE AnyWaveMontage>\n<Montage>\n');
for i = 1:size(channels,1)
    labels = strsplit(channels{i},'-');
    fprintf(fid, '\t<Channel name="%s">\n\t\t<type>SEEG</type>\n\t\t<reference>%s</reference>\n\t\t<color>black</color>\n\t</Channel>\n',...
            labels{1},labels{2});
%end
end
fprintf(fid, '</Montage>');
fclose(fid)
end
