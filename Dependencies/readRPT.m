function datamatrix = readRPT(filename)
%FB Version 08032021
%Read a Hanna Instruments .RPT file with titration data and extract an Nx2
%matrix with cols for volume, resulting pH. 
%INPUT: 
%   filename: valid file path for .RPT file
%OUTPUT:
%   datamatrix: Nx2 matrix of cumulative volume added (mL) and resulting pH
%NOTE: 
%   NO ERROR CHECKING is done for the file read
%------------------------------------------------------------------------
datamatrix = zeros(1000,2);                     %set empty data matrix
fid = fopen(filename, 'rt');                    %open file
readmatrix = 0;                                 %boolean for matrix data
indx = 1;
if fid ~= -1                                    %check for open file
   fline = fgetl(fid);                          %get 1st line
   while ~feof(fid)                             %while in file
      if contains(fline,'Volume[mL]')           %check for matrix header
         readmatrix = 1;                        %set matrix boolean true
      else                                      %not header line
         if readmatrix == 1                     %if in matrix
            if isempty(fline)                   %reached end of matrix
               break;                           %break out
            end
            C = regexp(fline,'\d+\.?\d*','match'); %Regex for floats
            vol = str2double(C{2});             %get float vol value
            pH = str2double(C{4});              %get float pH value
            datamatrix(indx,1) = vol;           %store in matrix (col 1)
            datamatrix(indx,2) = pH;            %store in matrix (col 2)
            indx = indx + 1;                    %advance matrix index
         end
      end
      fline = fgetl(fid);                       %read next line
   end
indices = (datamatrix(:,2) == 0);               %all rows with 0 in col 2
datamatrix(indices,:) = [];                     %delete those rows
end %end of file ID if statement
end

