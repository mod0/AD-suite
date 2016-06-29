function data = readFortranOutput(filename, row, col)
% To write the file in a format readable by this function, follow the
% example below:
% open(unit = 1, file = 'Pc1', form = 'unformatted', access = 'stream')
% write(unit=1) Pc(1, :)
% close(unit=1)
%
% Open the binary file in read mode
f = fopen(filename,'r');
% Read the data into an array, give the size in [row, col] format and 
% use double precision.
data = fread(f, [row,col], 'double');
end

