function newABmat = CombineABs( ABmat, tol )
%FB Version 08032021 
%PURPOSE: to combine buffers with similar pK values, prevent reduncancy in
% the curve fitting algorithm
%INPUT: 
%   ABmat: Nx2 matrix of acid conc (M) and pK value
%   tol: Scalar, similarity of pK, if two pK values in ABmat are less that
%    tol difference they are combined and concentrations added.
%OUTPUT: 
%   newABmat: Nx2 matrix acid conc (M) and pK with reduced pK values
%-------------------------------------------------------------------------

[row,~] = size(ABmat);                    %get number of AB pairs
newABmat = ABmat(1,:);               %new AB mat with 1st val
index = 1;                                %index for new ABmat
if row > 1
   for i=2:row
      if abs(ABmat(i,2) - ABmat(i-1,2)) < tol  %if next with tol of prev
         temp_conc = ABmat(i,1) + newABmat(index,1); %get new conc 
         newABmat(index,1) = temp_conc;      %set new value to new AB mat
      else
         index = index + 1;                  %advance index
         newABmat(index,:) = ABmat(i,:);     %pK values different, so use
      end 
   end
end

end

