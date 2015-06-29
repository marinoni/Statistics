function [H,PValue,Statistic,CriticalValue] = MCstatistics(x,varargin)

%MCstatistics performs godness-of-fit test for a number of distributions
%with a number of different methods on a one dimnesional data vector given in input.
%This routine does not work for matrices.
%
%The routine computes p-values and critical values basd on Monte-Carlo simulations 
%computed on the chosen method and significance value
%
%By calling the routine with all the parameters+1, not important the value, will display 
%debugging information 
%
%INPUTS:
%
% - X     : Data in input [Default, none]
% - ALPHA : Significance level [Default: 95%, i.e. 0.05]
% - DISTR : Nnull distribution [Default, 1]
%           * 1 - Normal
%           * 2 - Chi-squared
%           * 3 - T-Student
%           * 4 - Exponential
%           * 5 - Extreme Value
% - PARAM : Free parameters of the null distribution 
%           * 1 - Mean and standard deviation    [Default, mean and std of the data]
%           * 2 - Degrees of freedom             [Default, length of the data minus 1]
%           * 3 - Degrees of freedom             [Default, length of the data minus 1]
%           * 4 - Mean                           [Default, mean of the data]
%           * 5 - Location and Scale paramenters [Default, location and scale parameters of the data]
% - METHOD: Godness-of-fit implemented [Default, 1]
%           * 1 - Anderson-Darling
%           * 2 - Cramer-von Mises
%           * 3 - Kolmogorov-Smirnov     
% - MCTOL : Tolerance of the Monte-Carlo simulation on the p-Value [Default, 0.01]
%
%OUTPUTS:
%
% - H             : Result of the test (0 = the null hypothesys is accepted, 1 is rejected)
% - PVALUE        : P-Value of the statistics. Small values of PVALUE cast doubt on the validity of the null hypothesis
% - STATISTIC     : Statistic of the chosen method
% - CRITICALVALUE : Critical value of the chosen METHOD. When STATISTIC > CRITICALVALUE, the null hypothesis 
%                   can be rejected at a significance level of ALPHA.
%
%USAGE:
%
% [H,PVALUE,STATISTIC,CRITICALVALUE] = MCstatistics(X,ALPHA,DISTR,PARAM,METHOD,MCTOL)
%
%A. Marinoni, 15/05/2013

if ~isvector(x) || ~isreal(x)
    error('stats:lillietest:BadData',...
          'Input sample X must be a vector of real values.');
end

% Remove missing observations indicated by NaN's, and ensure that
% at least 4 valid observations remain.
x=x(~isnan(x));       % Remove missing observations indicated by NaN's.
x=x(:);               % Ensure a column vector.
n=length(x);
if n < 4
   error('stats:lillietest:NotEnoughData',...
         'Sample vector X must have at least 4 valid observations.');
end

defaults={0.05,1,0,1,0.01};
variab={'Alpha','Distr','Param','Method','MCtol'};
distr_poss={'Normal','Chi-squared','T-Student','Exponential','Extreme Value'};
param_poss=[2,1,1,1,2];
par_constraint=[0,0,2,-inf,0];
method_poss={'Anderson-Darling','Cramer Von Mises','Kolmogorov-Smirnov'};

lvarar=length(varargin)+1;
debug=0;
if lvarar>length(variab)+1
   lvarar=length(variab)+1;
   debug=1;
end

for i=1:lvarar-1
   if or(isempty(varargin{i}),~isnumeric(varargin{i}))
      eval(strcat([variab{i},'=',num2str(defaults{i}),';']));
   else
      eval(strcat([variab{i},'=[',num2str(varargin{i}),'];']));
   end
end
for i=lvarar:length(variab)
   eval(strcat([variab{i},'=',num2str(defaults{i})],';'));
end

if ~isscalar(Alpha) || ~(0<Alpha && Alpha<1)
   error('stats:MCstatistics:BadAlpha',...
         'Significance level Alpha must be a scalar between 0 and 1, abort');
end

indd=find(ismember([1:5],Distr));
if isempty(indd)
   error('stats:MCstatistics:BadDistribution: The required distribution is unknown, abort.');
end

indm=find(ismember([1:3],Method));
if isempty(indm)
   error('stats:MCstatistics:BadMethod: The required method is unknown, abort.');
end

% Calculate S(x), the sample CDF.
[sampleCDF,xCDF,n,emsg,eid] = cdfcalc(x);
if ~isempty(eid)
   error(sprintf('stats:MCstatistics:%s',eid),emsg);
end

%Compute null CDFs
switch Distr
   case 1 %normal
      if or(Param==0,param_poss(indd) ~= length(Param))
	 disp('Forcing mean and standard deviation of the null distribution to be those of the data')
	 Param=[mean(x) std(x)];
      end
      nullCDF = normcdf(xCDF, Param(1),Param(2));
   case 2 %chi-squared
      if or(Param==0,param_poss(indd) ~= length(Param))
         disp('Forcing degrees of freedom of the null distribution to be equal to the length of data minus one')
	 Param=n-1;   
      end
      nullCDF = chi2cdf(xCDF,Param);
   case 3 %T-Student
      if or(Param==0,param_poss(indd) ~= length(Param))
         disp('Forcing degrees of freedom of the null distribution to be equal to the length of data minus one')
	 Param=n-1;   
      end
      nullCDF = tcdf(xCDF,Param);
   case 4 %Exponential
      if or(Param==0,param_poss(indd) ~= length(Param))
         disp('Forcing degrees of freedom of the null distribution to be equal to the length of data minus one')
	 Param=mean(x);   
      end
      nullCDF = expcdf(xCDF, Param);
   case 5 %
      if or(Param==0,param_poss(indd) ~= length(Param))
         disp('Forcing locaion and scaleparameters of the null distribution to be given by type I maximum likelihood fit')
	 Param = evfit(x);
      end
      nullCDF = evcdf(xCDF, Param(1),Param(2));
   otherwise
       error(strcat(['stats:MCstatistics:BadDistribution. The null distribution must be one of ',distr_poss]));
end
   
%Check parameters values
if Param(end)<=par_constraint(Distr)
   strindx={'1st','2nd'};
   error(strcat(['stats:MCstatistics:BadParameter. The ',strindx{length(Param)},' free parameter of the ',distr_poss{Distr},...
   ' distribution must be larger than ',num2str(par_constraint(Distr))]))
end      
   
if ~isscalar(MCtol) || MCtol<=0
   error('stats:lillietest:BadMCReps',...
         'Monte-Carlo standard error tolerance MCTOL must be a positive scalar value.');
end

if debug
   disp('Summary')
   disp('-------')
   disp(strcat(['Method: ',method_poss{Method}]))
   disp(strcat(['Significance level: ',num2str(Alpha)]))
   disp(strcat(['Null distribution: ',distr_poss{Distr}]))
   disp(strcat(['Parameter(s): ',num2str(Param)]))
   disp(strcat(['Monte-Carlo tolerance on p-value: ',num2str(MCtol)]))
end
 
%Compute the sample statistic
[Statistic] = EDFdist(Method,sampleCDF,nullCDF,n);

% Compute the critical value and p-value on the fly using Monte-Carlo simulation.
[CriticalValue,PValue] = MCcrit(Statistic,n,Alpha,Distr,Param,Method,MCtol);

% Returning "H = 0" implies that we "Do not reject the null hypothesis at the
% significance level of alpha" and "H = 1" implies that we "Reject the null
% hypothesis at significance level of alpha."
H=double(Statistic > CriticalValue);

% ----------------------------------------------------------------------
function [Statistic] = EDFdist(Method,sampleCDF,nullCDF,n)
%
%This function takes in input the sample and null CDF as well as the method
%used to compute the deviation of the sample CDF from the null CDF, and returns
%the desired deviation.
%Allowed methods are:
%-1: Anderson-Darling
%-2: Cramer-Von Mises
%-3: Kolmogorov-Smirnov
%
%USAGE:
%
%[Statistic] = EDFdist(Method,sampleCDF,nullCDF,n)
%
%A. Marinoni, 15/05/2013

if ~ismember([1:3],Method)
   error('Statistics:BadMethod. Invalid method to compute the deviation from sample and null CDF: ',... 
          'choose 1, 2 or 3. Abort')
end   

what=0;
if what
   ind=find(nullCDF==0);
   nullCDF(ind)=eps*10.^ind;
   ind=find(nullCDF==1);
   nullCDF(ind)=1-eps*10.^(length(nullCDF)+1-ind);
end

delta1=sampleCDF(1:end-1)-nullCDF; % Vertical difference at jumps approaching from the LEFT.
delta2=sampleCDF(2:end)-nullCDF; % Vertical difference at jumps approaching from the RIGHT.
i=[1:n]';

switch Method

   case 1
   % Compute the test statistic of interest: T = ||S(x) - nullCDF(x)||_{L^2,weighed}.	
      Statistic=-n-sum(((2*i-1)/n).*(log(nullCDF)+log(1-nullCDF(n+1-i))));
      
      %The direct computation below works decently only when n is larger than a few hundred points
      %kk=linspace(1/n,n/(n+1),n);
      %kk=kk(:);
      %deltaCDF=(kk-nullCDF).^2;
      %deltaCDF=deltaCDF./nullCDF./(1-nullCDF);
      %deltaCDF=[0;deltaCDF;0];
      %nullCDF=[0;nullCDF;1];
      %Statistic=n*NC_integrate(deltaCDF,nullCDF,4,'quiet');

   case 2
   % Compute the test statistic of interest: T = ||S(x) - nullCDF(x)||_{L^2}.
   
      Statistic=1/12/n+sum(((2*i-1)/(2*n)-nullCDF).^2);
      %The direct computation below works decently only when n is larger than a few hundred points
      %kk=linspace(1/n,n/(n+1),n);
      %kk=kk(:);
      %deltaCDF=(kk-nullCDF).^2;
      %Statistic=NC_integrate(deltaCDF,nullCDF,1,'quiet');

   case 3
   % Compute the test statistic of interest: T = max|S(x) - nullCDF(x)|.
      deltaCDF=abs([delta1 ; delta2]);
      Statistic=max(deltaCDF); 
   
end   


% ----------------------------------------------------------------------
function [crit,p] = MCcrit(SampleStat,n,Alpha,Distr,Param,Method,MCtol)
%LILLIEMC Simulated critical values and p-values for Lilliefors' test.
%   [CRIT,P] = LILLIEMC(KSSTAT,N,ALPHA,DISTR,MCTOL) returns the critical value
%   CRIT and p-value P for Lilliefors' test of the null hypothesis that data
%   were drawn from a distribution in the family DISTR, for a sample size N
%   and confidence level 100*(1-ALPHA)%.  P is the p-value for the observed
%   value KSSTAT of the Kolmogorov-Smirnov statistic.  DISTR is 'norm', 'exp',
%   'or 'ev'. ALPHA is a scalar or vector.  LILLIEMC uses Monte-Carlo
%   simulation to approximate CRIT and P, and chooses the number of MC
%   replications, MCREPS, large enough to make the standard error for P,
%   SQRT(P*(1-P)/MCREPS), less than MCTOL.

vartol = MCtol^2;

crit = 0;
p = 0;
mcRepsTot = 0;
mcRepsMin = 1000;
yCDF = (0:n)'/n;
while true
    mcRepsOld = mcRepsTot;
    mcReps = ceil(mcRepsMin - mcRepsOld);
    statMC = zeros(mcReps,1);

    switch Distr

    % Simulate critical values for the normal   
    case 1
    
        for rep = 1:length(statMC)
            x = normrnd(Param(1),Param(2),n,1);
            xCDF = sort(x); % unique values, no need for ECDF
            nullCDF = normcdf(xCDF, mean(x), std(x)); % MLE fit to the data
            statMC(rep)=EDFdist(Method,yCDF,nullCDF,n);
        end
    
    % Simulate critical values for the chi2
    case 2
        
        for rep = 1:length(statMC)
            x = chi2rnd(Param,n,1);
            xCDF = sort(x); % unique values, no need for ECDF
            nullCDF = chi2cdf(xCDF,(mean(x)+var(x)/2)/2); % MLE fit to the data, dof_LeastSquares=(mean(x)+var(x)/2)/2 
            statMC(rep)=EDFdist(Method,yCDF,nullCDF,n);
        end
	
    % Simulate critical values for the t-stud
    case 3
        for rep = 1:length(statMC)
            x = trnd(Param,n,1);
            xCDF = sort(x); % unique values, no need for ECDF
            nullCDF = tcdf(xCDF,2*var(x)/(var(x)-1)); % MLE fit to the data, var(x)=dof/(dof-2)
            statMC(rep)=EDFdist(Method,yCDF,nullCDF,n);
        end

    % Simulate critical values for the exponential
    case 4

        for rep = 1:length(statMC)
            x = exprnd(Param,n,1);
            xCDF = sort(x); % unique values, no need for ECDF
            nullCDF = expcdf(xCDF, mean(x)); % MLE fit to the data
            statMC(rep)=EDFdist(Method,yCDF,nullCDF,n);
        end

    % Simulate critical values for the extreme value
    case 5
        
        for rep = 1:length(statMC)
            x = evrnd(Param(1),Param(2),n,1);
            xCDF = sort(x); % unique values, no need for ECDF
            phat = evfit(x); % MLE fit to the data
            nullCDF = evcdf(xCDF, phat(1), phat(2));
            statMC(rep)=EDFdist(Method,yCDF,nullCDF,n);
        end
    end

    critMC = prctile(statMC,100*(1-Alpha));
    pMC = sum(statMC > SampleStat) ./ mcReps;

    mcRepsTot = mcRepsOld + mcReps;
    crit = (mcRepsOld*crit + mcReps*critMC) / mcRepsTot;
    p = (mcRepsOld*p + mcReps*pMC) / mcRepsTot;

    % Compute a std err for p, with lower bound (1/N)*(1-1/N)/N when p==0.
    sepsq = max(p*(1-p)/mcRepsTot, 1/mcRepsTot^2);
    if sepsq < vartol
        break
    end

    % Based on the current estimate, find the number of trials needed to
    % make the MC std err less than the specified tolerance.
    mcRepsMin = 1.2 * (mcRepsTot*sepsq)/vartol;
end
