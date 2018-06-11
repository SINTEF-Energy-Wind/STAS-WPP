function Q = buildQij (Neb,nf,df,Vinf,I,Lu,Omega,r)
%
% We require the correlation functions between pairs of uz and ut
% turbulence components at each blade element. These are to be
% stored in the matrix Qij.  The columns (row index) are the
% functions of time offset tau (as this allows the entire matrix
% to be converted in one call to fft) and the columns are ordered
% (column index) as...
% Q[uz(el1  ,bl1),uz(el1  ,bl1),tau]
% Q[uz(el1  ,bl1),ut(el1  ,bl1),tau]
% Q[ut(el1  ,bl1),uz(el1  ,bl1),tau]
% Q[ut(el1  ,bl1),ut(el1  ,bl1),tau]
% Q[uz(el1  ,bl1),uz(el2  ,bl1),tau]
%               ...
% Q[ut(el1  ,bl1),ut(elNeb,bl1),tau]
% Q[uz(el1  ,bl1),uz(el1  ,bl2),tau]
%               ...
% Q[ut(el1  ,bl1),ut(elNeb,bl3),tau]
% Q[uz(el2  ,bl1),uz(el1  ,bl1),tau]
%               ...
% Q[ut(elNeb,bl1),ut(elNeb,bl3),tau]
%
% The total size is 4*nf-by-4*3*Neb^2.  Note that this takes
% advantage of the isotropy of the turbulence, such that we need
% consider only one reference blade, the relationships being
% identical for all three blades.
%
% (A thought.  For variations in speed, could this not be accounted
% for approximately by applying some representative variation, in 
% the time domain, to Omega as I compute Omega*tau?  Then I would 
% get a "smeared" spectrum like I would expect from time-domain
% simulations.)
%
% Version:        Changes:
% --------        -------------
% 14.02.2015      Adapted from the rotSpectrum subroutine in the
%                 modalAnalysis.f90 Fortran code.
% 02.08.2017      Adapted for complex step derivatives.
%
% Version:        Verification:
% --------        -------------
% 14.02.2015      Memo AN 15.12.19
% 02.08.2017      
%
% Inputs:
% -------
% Neb             : Number of elements per blade.
% nf              : Number of frequencies for which an accurate
%                   result is desired.  One-quarter of the number
%                   used in the computation.
% df              : Frequency bin width.
% Vinf            : Remote incoming windspeed.
% I               : Turbulence intensity = sig_u/Vinf.
% Lu              : Turbulence length scale.
% Omega           : Rotational speed.
% r               : Radius of each blade element.
%
% Outputs:
% --------
% Q               : Correlation matrix.

N = 4*nf; % Number of computed frequencies/tau intervals in order
          % to get nf useful, non-aliased intervals.
NN = 3*Neb;

Q = zeros(N,4*Neb*NN);

sig = I*Vinf;
third = 1/3;
twosiggam = 2*(sig^2)/gamma(third);
oneLu = 1.34*Lu;

T0 = 1/df;
dt = T0/N;
tau = [[0:dt:T0/2].';[-(T0/2-dt):dt:-dt].'];
Wt = Omega*tau;
twopi3 = 2*pi*third;
fourpi3 = 4*pi*third;
c2 = cos(twopi3);
s2 = sin(twopi3);
c4 = cos(fourpi3);
s4 = sin(fourpi3);
cwt = zeros(N,3);
swt = zeros(N,3);
cwt(:,1) = cos(Wt);
swt(:,1) = sin(Wt);
cwt(:,2) = cwt(:,1)*c2 - swt(:,1)*s2;
swt(:,2) = swt(:,1)*c2 + cwt(:,1)*s2;
cwt(:,3) = cwt(:,1)*c4 - swt(:,1)*s4;
swt(:,3) = swt(:,1)*c4 + cwt(:,1)*s4;

sz = -Vinf*tau;

for i1 = 1:Neb
   ib1 = 1;
   ival = 12*Neb*(i1-1);
   for i2 = 1:NN
      icol = ival + 4*(i2-1);
      ib2 = 1 + (i2>Neb) + (i2>2*Neb);

      if (ib1 == ib2)
         ind = 1;
      elseif (ib2 - ib1 == 1) || (ib2 - ib1 == -2)
         ind = 2;
      else
         ind = 3;
      end

      % Compute the relevant coordinates.
      sx = r(i2)*cwt(:,ind) - r(i1);
      if abs(real(r(i1)) - real(r(i2))) < ...
                      0.001*(real(r(i1)) + real(r(i2)))
         % Prevent later divide-by-zero error when tau = 0.
         ii = abs(real(tau)) < 0.5*real(dt);
         sx(ii) = sx(ii) - 0.005*r(i1);
      end
      sy = r(i2)*swt(:,ind);
      s = sqrt(sx.^2 + sy.^2 + sz.^2);
      sLu = s/oneLu;

      Qss = twosiggam*((0.5*sLu).^third).*besselk(third,sLu);
      dQss = -(twosiggam/oneLu)*((0.5*sLu).^third).*besselk(-2*third,sLu);
      dQsss = dQss./s;

      % Intermediate correlations.
      Qyy = Qss + 0.5*s.*dQss - 0.5*(sy.^2).*dQsss;
      Qyx = -0.5*sy.*sx.*dQsss;
      Qzx = -0.5*sz.*sx.*dQsss;
      Qzy = -0.5*sy.*sz.*dQsss;
      % Qyz = Qzy;

      % Desired correlations.
      Q(:,icol+1) = Qss + 0.5*s.*dQss - 0.5*(sz.^2).*dQsss;
      Q(:,icol+2) = -swt(:,ind).*Qzx + cwt(:,ind).*Qzy;
      Q(:,icol+3) = Qzy; % Qyz;
      Q(:,icol+4) = -swt(:,ind).*Qyx + cwt(:,ind).*Qyy;

   end
end


