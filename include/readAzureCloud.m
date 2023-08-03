function [matrix, info] = readAzureCloud(fileName)
fid = fopen(fileName);
info = fread(fid, 2, 'int32');
B = fread(fid, info(1)*info(2) * 3, 'int16',0,'l');
matrix = reshape(B, 3, length(B)/3)';
fclose(fid);
end
