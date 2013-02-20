function  rawData = ITKRawFileReader(filename)

% Reads a raw file
% Inputs: filename, image size, image type, byte ordering, skip
% If you are using an offset to access another slice of the volume image
% be sure to multiply your skip value by the number of bytes of the 
% type (ie. float is 4 bytes).
% Inputs: filename, image size, pixel type, endian, number of values to
% skip.
% Output: image
%example: rawData = ITKRawFileReader('bin\flipped.raw',[800 600
%40],'float32','l',0);
% Author : Amit Mukherjee


[imSize, type] = readMetaHeader(filename);

if numel(imSize) ~= 3
    error('unable to read');
end

endian = 'l';
skip = 0;
fid = fopen(filename,'rb',endian);
if (fid < 0)
    fprintf('Filename %s does not exist\n',filename);
    rawData = -1;
else
    status = fseek(fid,skip,'bof');
	if status == 0
        rawData = fread(fid,prod(imSize),type);
        fclose(fid);
        if (length(imSize) == 3)
            slices = length(rawData)/imSize(1)/imSize(2);
            imSize(3) = slices;
        end
        rawData = reshape(rawData,imSize);
	else
        rawData = status;
	end
end
rawData = permute(rawData,[2,1,3]);

function [imsize, type] = readMetaHeader(filename)
% [name, ext] = strtok(filename,'.');
name = filename(1:end-4);
fid = fopen([name,'.mhd'],'r');
if fid < 0 
    disp(['Header file : ', [name,'.mhd'], ' not found']);
    imsize = [];
    type = [];
    return;
end

progress = 0;
while 1
    tline = fgetl(fid);
    if ~ischar(tline),   break,   end
    disp(tline);
    [a, b] = strtok(tline,'=');
    if (~isempty(b))
        if strcmpi(a,'DimSize ')
           imsize = str2num(b(2:end)); 
           progress = progress + 1;
        elseif strcmpi(a,'ElementType ')
            if strcmp(b(3:end), 'MET_FLOAT')
                type = 'float32';
                progress = progress + 1;
            end
        end
    end
end
fclose(fid);
if progress ~= 2
    disp(['Header file : ', [name,'.mhd'], ' Parameter error']);
    disp(imsize);
    disp(type);
    imsize = [];
    type = [];
end
