function f = getFrequencies (Nf,df)
%
% If uniform frequency bins are used, then the waves become
% periodic at 1/df.  Smarter to sample randomly or with equal
% areas to avoid this.  
%
% Here we employ a simplified version of Shinozuka's method,
% assuming a sine function for the CDF.
%

fmin = df;
fmax = Nf*df;
f1 = fmin + (fmax-fmin)*(asin(1.6*rand(Nf,1) - 0.8) ...
   / (asin(0.8) - asin(-0.8)) + 0.5);

[f,j] = sort(f1,'ascend');


