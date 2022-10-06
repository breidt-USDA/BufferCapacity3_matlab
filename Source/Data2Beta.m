function res = Data2Beta( acidTitration, baseTitration, molarityHCl, ...
   molarityNaOH, initialvol, minDeltapH)
%Version FB_08242017
%DESCRIPTION
%Takes acid (HCl) and base (NaOH) titration and generates buffer capacity
%  (BC) curve.
%PARAMETERS
%  acidTitration = N1x2, volume added (cumulative) in ml, resulting pH
%  baseTitration = N2x2, volume added (cumulative) in ml, resulting pH
%  molarityHCl = concentration of acid (HCl) in M/L used for titration
%  molarityNaOH = concentration of base (NaOH) in M/L used for titration
%  vol = initial volume of solution to be titrated (L)
%  minDeltapH = minimum change in pH for an added aliquote of acid or base
%RETURNS
%  res = (Nx2) 
%   Col1: pH value after at each titration data point
%   Col2: dn/dpH = change in equivalents (acid or base) Molar change  
%     divided by the change in pH, and corrected for volume.
%OUTPUT CSV FILES (supressed in this version)
%  acid_data.csv, contains data on acid titrtion
%  base_data.csv, contains data on base titration
%  BC_curve.csv, contains data for BC curve, including indices for each
%     data point (this shows which points removed). 
%DEPENDENCIES
%  Matlab bulit-in Array2Table  
%NOTE: First 2 junction data points are removed because delta pH change is
%  exaggerated. 
%NOTE: derivative values resulting in little or no pH change on addtion
%  of acid (pH not going down) or base (pH not going up) are removed 
%NOTE: BC is defined as the change in the MOLARITY of acid/base resulting
%  from additon of acid/base divided by the resulting change in pH 
%------------------------------------------------------------------------

%report initial pH values for acid titraiton and base titration
res.init_titration_vol = initialvol;   %initial titration vol in L
res.init_pH_acid = acidTitration(1,2); %initial pH, acid titraiton
res.init_pH_base = baseTitration(1,2); %initial pH, base titration

%setting for minimum delta pH for each titration step
removeVal = minDeltapH;                %removeal criterion for delta pH
 
%get BC results for each titration and add index values for each data set
res.acid = tcurve(acidTitration, molarityHCl, initialvol);
res.acid(:,7) = (length(res.acid):-1:1)'; %num backwards for acid data
res.base = tcurve(baseTitration, molarityNaOH, initialvol); 
res.base(:,7) = (1:length(res.base))+length(res.acid); %num from acid val

%convert units for acid/base and format output data 
res.acid(:,1) = round((res.acid(:,1)*1000),4); %L to ml total vol
res.acid(:,2) = round(res.acid(:,2),4);  %pH
res.acid(:,3) = round(res.acid(:,3),4);  %delta pH
res.acid(:,4) = round((res.acid(:,4) * 1000000),0);  %L to ul vol added
res.acid(:,5) = round((res.acid(:,5) * 1000),4);  %M/L to mM/L delta acid 
res.acid(:,6) = round(res.acid(:,6),4);  %BC value
res.acid(:,7) = round(res.acid(:,7),0);  %Index
res.base(:,1) = round((res.base(:,1)*1000),4); %L to ml total vol
res.base(:,2) = round(res.base(:,2),4);  %pH
res.base(:,3) = round(res.base(:,3),4);  %delta pH
res.base(:,4) = round((res.base(:,4) * 1000000),0);  %L to ul vol added
res.base(:,5) = round((res.base(:,5) * 1000),4);  %M/L to mM/L delta base 
res.base(:,6) = round(res.base(:,6),4);  %BC value
res.base(:,7) = round(res.base(:,7),0);  %Index
   
%print output acid/base csv files 
% Tab = array2table(res.acid, 'VariableNames', {'Total_Vol_ml','pH', ...
%    'Delta_pH','Vol_Added_ul','Delta_Acid_mMperL','BC_value','Index'});
% writetable(Tab, 'acid_data.csv');
% Tab = array2table(res.base, 'VariableNames', {'Total_Vol_ml','pH', ...
%    'Delta_pH','Vol_Added_ul','Delta_Base_mMperL','BC_value','Index'});
% writetable(Tab, 'base_data.csv');

%remove 1st two pts (start of titration 0 vals + 1)
temp.acid = res.acid;         %copy acid data to adjust for BC curve
temp.acid(1:2,:) = [];        %remove junction values based on 0 vol
temp.base = res.base;         %copy base data to adjust for BC curve
temp.base(1:2,:) = [];        %remove junction values based on 0 vol

%flip the acid titration and append base titration
BC = [flipud(temp.acid);temp.base];

%use logical indexing to remove rows w/specific condition (col 3 < val)
TF = BC(:,3) < removeVal;     %generate True-False vector (TF)
BC(TF,:) = [];                %adjust matrix by logical val

%reformat BC data
BC(:,1) = BC(:,2);            %pH values to 1st col
BC(:,2) = BC(:,6);            %Calc BC vals to 2nd col
BC(:,3) = BC(:,7);            %index values to 3rd col
BC(:,4:7) = [];               %delete rest of matrix

%output BC curve csv file
% Tab = array2table(BC,'VariableNames', {'pH','BC','Index'});
% writetable(Tab, 'BC_curve.csv');
BC(:,3) = [];                 %delete index column
res.BCcurve = BC;             %output for BC curve (Nx2): BC, pH

%===============Supporting_Functions======================================
function rval = tcurve(titration, molarity, initvol) 
   %titration = Nx2: ml added (accum), resulting pH
   %molarity = molarity of titrant
   %initial vol in liters
   %----------------------------------------------------------------------
   vols = titration(:,1)/1000;      %convert accum added vols to liters(L)
   nvals = length(vols);            %get number of data points
   rval = zeros(nvals,6);           %set up return matrix
   rval(:,1) = vols + initvol;      %convert to total accum vol (L)
   rval(:,2) = titration(:,2);      %assign pH values for each step
   for i=2:nvals                    %for all additions of titrant:
      rval(i,3) = abs(rval(i,2) - rval(i-1,2)); %delta pH 
      rval(i,4) = vols(i) - vols(i-1); %delta vol of added titrant
      rval(i,5) = rval(i,4)*molarity/rval(i,1); %delta M/L titrant in sol
      if rval(i,3) ~= 0             %check for divide by zero error
         rval(i,6) = rval(i,5)/rval(i,3); %calc BC = dC/dpH
      end
   end
end 


end


