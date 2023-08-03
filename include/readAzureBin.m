function [matrix, info] = readAzureBin(fileName)
fid = fopen(fileName);
info = fread(fid, 2, 'int32');
B = fread(fid, info(1)*info(2), 'uint16',0,'l');
matrix = reshape(B,info(1),info(2))';
fclose(fid);
end
