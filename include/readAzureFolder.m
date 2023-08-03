function frames = readAzureFolder(folderName, raw)
%To plot azure CLoud use:
%Cloud data is N by 3
%plot(cloud(:,1), cloud(:,3), -cloud(:,2))

frames = [];
fid = fopen([folderName, '/log.txt']);

while(~feof(fid))
    row = textscan(fid, '%d16 %s %d16 %s %d16 %s\r', 'Delimiter',',');
    if (length(row) == 6)
        if (~isempty(row{2}))
            frame.index = row{1};
            frame.ts = row{2}{1};
            frame.rgb.index = row{3};
            frame.rgb.ts = row{4}{1};
            frame.rgb.data = imread([folderName,'/color',num2str(frame.rgb.index),'.png']); 
            if (raw)
                frame.rgb.data = flipud(frame.rgb.data);
            end
            frame.depth.index = row{5};
            frame.depth.ts = row{6}{1};
            frame.depth.data = readAzureBin([folderName,'/depth',num2str(frame.depth.index),'.bin']);
%             frame.depth.data = readAzureBin([folderName,'/depth',num2str(frame.depth.index),'.bin']); 
            frame.cloud.data = readAzureCloud([folderName,'/cloud',num2str(frame.depth.index),'.bin']); 
            frames = [frames;frame];
        end
    end
end
fclose(fid);

end

