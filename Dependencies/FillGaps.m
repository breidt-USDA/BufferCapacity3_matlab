function newBC = FillGaps(BCcurve, mingap, increment)
%FB Version 08042021
%Takes a BC curve with gaps, generates a gap list and an equation for
%a line that extends across the gap. Next, using the increment value, new
%X values (pH vals) are spaced evenly at the increment values across the  
%gap, and Y values on the line are added for each X. 
%
%Input: 
%  BCcurve = Nx2, with pH and BC value
%  mingap = gap size above which a fill is needed
%  increment = pH increment 
%Output:
%  newBC = new Nx2 BC curve (pH, BC). 
%
%NOTE: if no gaps, return original curve
%-------------------------------------------------------------------------

%scroll through BC curve and find all gaps larger than mingap
gaplist = FindGaps(BCcurve,mingap);    %identify the gaps
if all(gaplist(:) == 0)                %check to see if gaps existed
   newBC = BCcurve;                    %just return original function
   return;                             %end function if no gaps
end

%setup for new gap curve
[rows,~] = size(gaplist);              %get number of gaps in the list

newBC = BCcurve;                       %set up new BC curve, return matrix

%process gap list and add fills to newBC curve
for i=1:rows                           %for each gap
   tempStart = gaplist(i,2);           %get the starting pH
   tempEnd = gaplist(i,3);             %get the end pH 
   %generate vector with x values for gap
   xvalvec = GetXvalVec(tempStart,tempEnd,increment);
   Xvals = [gaplist(i,2) gaplist(i,3)]; %pH vals for gap start and end
   Yvals = [gaplist(i,4) gaplist(i,5)]; %BC vals for gap start and end
   %get the slop and intercept of the line between gap start-end points
   lineParams = GetLineParams(Xvals, Yvals);
   %now, get Y values for points on line at each xvalvec
   yvalvec = (lineParams(1) .* xvalvec) + lineParams(2);
   %make xy vectors into Nx2
   tempmat = [xvalvec' yvalvec']; 
   newBC = [newBC; tempmat];  
end

%sort new BC curve with filled gap data by pH values
newBC = sortrows(newBC,1);
end %end of FillGaps 

%SUBFUNCTIONS
function xvec = GetXvalVec(startval, endval, increment)
   %get a vector of x values (pH values) across gap for each increment
   tempval = startval;                 %starting pH
   indx = 0;                           %temp index
   while true                          %loop till end pH reached
      tempval = tempval + increment;   %increment pH
      if tempval < endval              %check for end value
         indx = indx + 1;              %next index for vector
         xvec(indx) = tempval;         %assign next pH value
      else
         break;                        %end pH reached, finish loop
      end
   end
end %end of GetXvalVec

function parameters = GetLineParams(Xs, Ys)
   %generate slope, intercept parameters for line from two XY points 
   slope = (Ys(2) - Ys(1))/(Xs(2) - Xs(1)); %calc. slope
   intercept = Ys(1) - (slope*Xs(1));  %calc intercept     
   parameters = [slope intercept];     %return line parameters
end %end of GetLineParams

function GapList = FindGaps(BCcurve, minGapSize)
%Function takes a BC curve, scans X values for gaps then makes a list of
%each gap bigger than the min gap size.
%Arguments: 
%  BCcurve, Nx2: cols of pH, BC val
%  minGapSize, scalor: smallest difference in pH that needs to be fixed 
%Ouptup: 
%  GapList, Nx5, cols:
%     1. index value of gap start (integer)
%     2. pH value at gap start index
%     3. pH value at gap end index (next index from gap start index
%     4. BC value at gap start index
%     5. BC value at gap end index
%NOTE: if no gaps return a 1 row gap list with all zeros
%-----------------------------------------------------------------------

%get starting vectors and values
Xvals = BCcurve(:,1);                     %X val vector
Yvals = BCcurve(:,2);                     %Y val vector
len = length(Xvals);                      %length (scalor)
GapVec = zeros(len-1,3);                  %Nx2 of start index and gaps

%deterime gap size
for i = 2:len                             %for each gap in list
   GapVec(i-1,1) = i-1;                   %start index for each gap
   GapVec(i-1,2) = i;                     %end index for each gap
   GapVec(i-1,3) = Xvals(i) - Xvals(i-1); %get each gap value
end

%sort gap vector largest gap to smallest
GapVec = sortrows(GapVec,-3);             %largest gaps at begining
ngaps = 0;
%build return matirx with only gaps that are larger than mingap size

for j = 1:len-1                           %scroll through gap list
   if GapVec(j,3) > minGapSize            %check gap size
      start_index = GapVec(j,1);          %if too big, start val index
      end_index = GapVec(j,2);            %assign end val index
      GapList(j,2:3) = [Xvals(start_index) Xvals(end_index)]; %Xvals
      GapList(j,4:5) = [Yvals(start_index) Yvals(end_index)]; %Yvals
      GapList(j,1) = start_index;         %start index for gap
      ngaps = ngaps + 1;
   else
      break;                              %remaining gaps too small
   end
end

if ngaps == 0                             %if no gaps
   GapList = zeros(1,5);                  %return one row of zeros     
end

end %end of FindGaps

