function indices = splitvec( byN, vecsize)
%FB Version 08032021
%Generate index values that are (approximately) evenly spaced across for
% a vector of a given size
%NOTE: always has 1 as a starting point
%INPUT: 
%   byN: number of indices needed (scalar)
%   vecsize: length of vector
%RETURNS:
%   indices: 1xN vector of index values
%NOTE: 
%   NO ERROR CHECKING of input values for byN within vecsize
%-------------------------------------------------------------------------

indices = zeros(1,byN);             %initialize return vector
increment = vecsize/(byN + 1);      %calc increment between values
previous = 0;                       %initialize previous index value
for i=1:byN                         %for each index in the vector
    indices(i) = round(previous + increment); %get int value for next index
    previous = indices(i);          %reset previous value
end

end