function res = SimplexBCPK_DF5( ObsXY, ABmat, minConc, ...
   pK_tol, LB, UB)
%FB version 01032022
%Optimize BC solution from observed BC data for a given set of pK values
%   at a given salt concentration. This version uses fmincon to fit 
%   concentration and pK values. Starting values are the BC values in 
%   ABmat, with target_pH being the mean initial pH form the titrations.
%NOTE that: 
% 1. Neither the ionic str. nor NaCl percent  were included in the beta 
%    model
%   calculations because the salt effect on pK is observed during
%   titration. 
%PARAMETERS
%   ObsXY: Nx2 matrix (pH, BCval)
%   ABmat: Nx2 matrix (Conc M/L, pH) 
%   minConc: min comcentration (M) for buffer (< val, buffer deleted)
%   pK_tol: minimum pH difference between buffers (<val, buffers combined)
%   LB: lower pH bound for optimization (scalar)
%   UB: upper pH bound for optimization (scalar)
%RETURNS
%   Matlab structure including:
%   res.Obs = ObsXY;                       %Nx2 col vec (pH, BC)
%   res.PredX = X values only 
%   res.Pred = PredData;                   %Nx2 vector (pH, BC)
%   res.ABmat = ABmat_pred;                %Nx2 output (pred conc, pKvals)
%   res.BCmat = Conc2Beta(res.ABmat);      %Nx2 output (pred beta, pKvals)
%   res.pH = CalcpH_AB(ABmat_pred,salt,0); %predicted pH
%   res.SSE = fval;                     %scalar (final sum sq err)
%   res.ExitFlag = exitflag;            %fminsearch flag (1,0,-1: scalar)
%   res.Evals = optvals.funcCount;      %number of function calls (scalar)
%DEPENDENCIES: 
% BCcurve = BetaModel_AB(Conc,BCmat,NaClpercent,pHvec)
%   Conc = col vector of buffer concentrations (M/L)
%   pKa = col vector of corresponding pKa values
%   NaClpercent = Percent NaCl in solution (needed for ionic str. calc)
%   pHvec = col vector of pH values at which to determine B val 
%   RETURNS: Buffer Capacity Curve (Nx2: pH, BC)  
% ABmat = CombineABs(ABmat, pK_tol); 
%   ABmat = Nx2, conc, pH
%   pK_tol = scalar (if difference in pK values is smaller, combine)
%   RETURNS: Nx2 ABmatrix that combines concentrations for similar pKs
% BCmat = Conc2Beta(ABmatrix)
%   ABmatrix = Nx2, conc, pK, or pH
%   RETURNS: Nx2 buffer capacity, pK matrix from Nx2 conc, pK matrix
%   NOTE: does not use %NaCl, only returns a BC value for pH value
%   DEPENDENCY: BetaModel_AB:
%NOTES:
% 1. minConc = 0.001, delete any buffer, pK pair with conc. below 1 mM
% 2. pK_tol = 0.1, combine buffer, pK pairs with pKs < tol apart. 
% 3. fmincon is constrained for each param, using limits for concentration 
%     as zero to 1 M/L, and each pK between min and max values for 
%     titration data. 
% 4. fmincon is using "interior_point" algorithm in this version
%*************************************************************************
ObsX = ObsXY(:,1);                     %pH values (col)
ObsY = ObsXY(:,2);                     %BC values (col)
zeroflag = false;                      %bool for check for no buffers
PredX = ObsX;                          %X value vector
salt = 0;                              %all calc done without IS correction
TF = abs(ABmat(:,1)) <= minConc;       %condition for minimum conc in AB
ABmat(TF,:) = [];                      %remove values with low conc
if isempty(ABmat)                      %check for empty matrix
   ABmat = [0 0];                      %if empty, just use zeros
   zeroflag = true;                    %set a flag so fit only water?
else
   ABmat = CombineABs(ABmat, pK_tol);  %combine rows w/similar pKs
end

%set up boundry values for fmincon constraints on each parameter value
[row,col] = size(ABmat);               %get size of ABmat
LowerBounds = zeros(row,col);          %LB matrix, conc values remain zero       
LowerBounds(:,2) = LB;                 %LB pKs have min of titration curve
UpperBounds = zeros(row,col);          %UB matrix
UpperBounds(:,2) = UB;                 %UB pKs have max of titration curve
UpperBounds(:,1) = 1;                  %UB conc values have max of 1 M/L 
lbvec = ABmat2paramvec(LowerBounds);   %convert matrix to vector (subfunc)
ubvec = ABmat2paramvec(UpperBounds);   %convert matrix to vector (subfunc)

%set up parameters for fmincon optimization
paramvec = ABmat2paramvec(ABmat);      %convert to Nx2 to 2N vector
fhan = @(paramvec)CalcSSE(paramvec, salt, ObsX, ObsY); %anon. function SSE
options = optimoptions('fmincon',...   %define options for optimization
   'Display','off',...
   'PlotFcn','optimplotfval');         %plot function value      

%optimize with constraints
[params, fval , exitflag, optvals]  = fmincon(fhan,paramvec,[],[],...
   [],[],lbvec,ubvec,[],options);

%regenerate AB matrix
ABmat_pred = paramvec2ABvals(params);  %convert linear vector to ABmat

%adjust AB matrix (remove very small or neg values and combine like pKs)
pK_tol = 0.3;                          %combine similar pKs 
TF = abs(ABmat_pred(:,1)) <= minConc;  %condition for minimum conc in AB
ABmat_pred(TF,:) = [];                 %remove values with low conc
ABmat_pred = CombineABs(ABmat_pred, pK_tol); %combine rows w/similar pKs 

%set pred data
if zeroflag
   ABmat_pred = [0 0];
   PredData = BetaModel_AB(ABmat_pred,salt, PredX); %returns Nx2 (pH,BC)
else
   %Predicted data with BetaModel
   PredData = BetaModel_AB(ABmat_pred,salt, PredX); %returns Nx2 (pH,BC)
end

%Record results
res.Obs = ObsXY;                       %Nx2 (pH vector, BC)
res.PredX = PredX;                     %Nx1 X vals (pH vector)
res.Pred = PredData;                   %Nx2 (pH vector, BC) 
res.ABmat = ABmat_pred;                %Nx2 output (pred conc, pKvals)
if sum(ABmat_pred(:,1) == 0 > 1)       %dont call with a zero conc row
   res.BCmat = [0 0];                  %return 1x2 zero matrix
else
   res.BCmat = Conc2Beta(ABmat_pred);  %get beta values from AB matrix
end
res.SSE = fval/20;                     %scalar (sum sq err, recorrect)
res.exitflag = exitflag;               %fminsearch flag (1,0,-1: scalar)
res.funcCount = optvals.funcCount;     %number of function calls (scalar)

%==SUBFUNCTIONS===========================================================
%Sum Squared Error term for predicted data with BetaModel 
function SSE = CalcSSE(paramvec, salt, ObsX, ObsY)
   %to help prevent rounding errors multiply by 20
   ABvalues = paramvec2ABvals(paramvec);  %convert linear vec to ABmat
   Pred = BetaModel_AB(ABvalues, salt, ObsX); %pred for given conc  
   SSE = 20*sum((ObsY-Pred(:,2)).^2);     %sum squared errors for BC vals
end

%Convert AB matrix to linear vector conc;pK in order
function paramvec = ABmat2paramvec( ABvalues)
   paramvec = [ABvalues(:,1);ABvalues(:,2)];  %conc first, then pKs
end

%Convert linear conc;pK params back to AB matrix 
function ABvals = paramvec2ABvals( paramvec ) 
   nvals = length(paramvec)/2;            %get index for 1/2 of vector
   ABvals(:,1) = paramvec(1:nvals);       %conc vals in 1st col
   ABvals(:,2) = paramvec(nvals+1:end);   %pK vals in 2nd col
end

end

