function Q = buildQijTower (Neb,Nut,nf,df,Vinf,I,Lu,Omega, ...
                            r,zut_g,sgr_g,Tg_r,Ty_g)
%
% High-wind load cases require that turbulent loads on the tower
% are included in the dynamic analysis.  This in turn requires
% that we augment the cross-spectral matrix of turbulence to 
% include points along the (stationary) tower.  We need then to
% formulate cross-correlations between uz,ut components on the
% three blades (the spanwise component being negligible), and
% ux,uy components on the tower (the uz component being 
% negligible).
%
% The resulting correlations are stored with columns ordered as
% Q[uz(el1  ,bl1),ux(el1  ,tow),tau]
% Q[uz(el1  ,bl1),uy(el1  ,tow),tau]
% Q[ut(el1  ,bl1),ux(el1  ,tow),tau]
% Q[ut(el1  ,bl1),uy(el1  ,tow),tau]
% Q[uz(el1  ,bl1),ux(el2  ,tow),tau]
%               ...
% Q[ut(el1  ,bl1),uy(elNut,tow),tau]
% Q[uz(el2  ,bl1),ux(el1  ,tow),tau]
%               ...
% Q[ut(elNeb,bl1),uy(elNut,tow),tau]
% Q[uz(el1  ,bl2),ux(el1  ,tow),tau]
%               ...
% Q[ut(elNeb,bl3),uy(elNut,tow),tau]
% Q[ux(el1  ,tow),uz(el1  ,bl1),tau]
%               ...
% Q[uy(elNut,tow),ut(elNeb,bl3),tau]
% Q[ux(el1  ,tow),ux(el1  ,tow),tau]
%               ...
% Q[uy(el1  ,tow),uy(el1  ,tow),tau]
% Q[ux(el1  ,tow),ux(el2  ,tow),tau]
%               ...
% Q[uy(el1  ,tow),uy(elNut,tow),tau]
% Q[ux(el2  ,tow),ux(el1  ,tow),tau]
%               ...
% Q[uy(elNut,tow),uy(elNut,tow),tau]
%
% The total size is 4*nf-by-4*(2*Nut*3*Neb + Nut^2).
%
% Version:        Changes:
% --------        -------------
% 05.10.2015      Original code.
% 04.08.2017      Adapted for complex step derivatives.
%
% Version:        Verification:
% --------        -------------
% 05.10.2015      Memo AN 15.12.68
% 04.08.2017      
%
% Inputs:
% -------
% Neb             : Number of elements per blade.
% Nut             : Number of turbulence points on the tower.
% nf              : Number of frequencies for which an accurate
%                   result is desired.  One-quarter of the number
%                   used in the computation.
% df              : Frequency bin width.
% Vinf            : Remote incoming windspeed.
% I               : Turbulence intensity = sig_u/Vinf.
% Lu              : Turbulence length scale.
% Omega           : Rotational speed.
% r               : Radius of each blade element.
% zut_g           : z coordinates of turbulence points on tower, in
%                   global coordinates.
% sgr_g           : Location of rotorplane (hub) center in global
%                   coordinates.
% Tg_r            : 3-by-3 transform between global and rotorplane.
% Ty_g            : 3-by-3 transform between yaw and global.
%
% Outputs:
% --------
% Q               : Correlation matrix.

N = 4*nf; % Number of computed frequencies/tau intervals in order
          % to get nf useful, non-aliased intervals.
NN = 3*Neb;

Q = zeros(N,4*(6*Nut*Neb + Nut^2));

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

sxy_g = Ty_g*[-Vinf*tau.';zeros(1,N);zeros(1,N)];

for i1 = 1:NN

   ival = 4*Nut*(i1-1);

   if (i1 == 1)
      ibl = 1;
      Tgh21 = -swt(1,ibl)*Tg_r(1,1) + cwt(1,ibl)*Tg_r(2,1);
      Tgh22 = -swt(1,ibl)*Tg_r(1,2) + cwt(1,ibl)*Tg_r(2,2);
      Tgh23 = -swt(1,ibl)*Tg_r(1,3) + cwt(1,ibl)*Tg_r(2,3);
      Tgh31 = Tg_r(3,1);
      Tgh32 = Tg_r(3,2);
      Tgh33 = Tg_r(3,3);
   elseif (i1 == Neb+1)
      ibl = 2;
      Tgh21 = -swt(1,ibl)*Tg_r(1,1) + cwt(1,ibl)*Tg_r(2,1);
      Tgh22 = -swt(1,ibl)*Tg_r(1,2) + cwt(1,ibl)*Tg_r(2,2);
      Tgh23 = -swt(1,ibl)*Tg_r(1,3) + cwt(1,ibl)*Tg_r(2,3);
      Tgh31 = Tg_r(3,1);
      Tgh32 = Tg_r(3,2);
      Tgh33 = Tg_r(3,3);
   elseif (i1 == 2*Neb+1)
      ibl = 3;
      Tgh21 = -swt(1,ibl)*Tg_r(1,1) + cwt(1,ibl)*Tg_r(2,1);
      Tgh22 = -swt(1,ibl)*Tg_r(1,2) + cwt(1,ibl)*Tg_r(2,2);
      Tgh23 = -swt(1,ibl)*Tg_r(1,3) + cwt(1,ibl)*Tg_r(2,3);
      Tgh31 = Tg_r(3,1);
      Tgh32 = Tg_r(3,2);
      Tgh33 = Tg_r(3,3);
   end

   s1_g = sgr_g + (Tg_r.')*[r(i1)*cwt(1,ibl);r(i1)*swt(1,ibl);0];

   for i2 = 1:Nut

      icol = ival + 4*(i2-1);

      s2_g = [0;0;zut_g(i2)] + sxy_g;
      s12 = (s2_g - s1_g).';
      s = sqrt(s12(:,1).^2 + s12(:,2).^2 + s12(:,3).^2);
      sLu = s/oneLu;

      Qss = twosiggam*((0.5*sLu).^third).*besselk(third,sLu);
      dQss = -(twosiggam/oneLu)*((0.5*sLu).^third).*besselk(-2*third,sLu);
      dQsss = dQss./s;

      % Intermediate correlations in global coordinates.
      Qxx = Qss + 0.5*s.*dQss - 0.5*(s12(:,1).^2).*dQsss;
      Qxy = -0.5*s12(:,1).*s12(:,2).*dQsss;
      Qyx = Qxy;
      Qyy = Qss + 0.5*s.*dQss - 0.5*(s12(:,2).^2).*dQsss;
      Qzx = -0.5*s12(:,3).*s12(:,1).*dQsss;
      Qzy = -0.5*s12(:,3).*s12(:,2).*dQsss;

      % Transform the blade components from the global to the
      % uz,ut (hub z,y) frame.
      Q(:,icol+1) = Tgh31*Qxx + Tgh32*Qyx + Tgh33*Qzx; % uz,ux
      Q(:,icol+2) = Tgh31*Qxy + Tgh32*Qyy + Tgh33*Qzy; % uz,uy
      Q(:,icol+3) = Tgh21*Qxx + Tgh22*Qyx + Tgh23*Qzx; % ut,ux
      Q(:,icol+4) = Tgh21*Qxy + Tgh22*Qyy + Tgh23*Qzy; % ut,uy

   end

end

ival0 = 4*Nut*NN;

for i1 = 1:Nut

   ival = ival0 + 4*NN*(i1-1);

   s1_g = [0;0;zut_g(i1)];

   for i2 = 1:NN

      icol = ival + 4*(i2-1);

      if (i2 == 1)
         ibl = 1;
         Tgh21 = -swt(:,ibl)*Tg_r(1,1) + cwt(:,ibl)*Tg_r(2,1);
         Tgh22 = -swt(:,ibl)*Tg_r(1,2) + cwt(:,ibl)*Tg_r(2,2);
         Tgh23 = -swt(:,ibl)*Tg_r(1,3) + cwt(:,ibl)*Tg_r(2,3);
         Tgh31 = Tg_r(3,1);
         Tgh32 = Tg_r(3,2);
         Tgh33 = Tg_r(3,3);
      elseif (i2 == Neb+1)
         ibl = 2;
         Tgh21 = -swt(:,ibl)*Tg_r(1,1) + cwt(:,ibl)*Tg_r(2,1);
         Tgh22 = -swt(:,ibl)*Tg_r(1,2) + cwt(:,ibl)*Tg_r(2,2);
         Tgh23 = -swt(:,ibl)*Tg_r(1,3) + cwt(:,ibl)*Tg_r(2,3);
         Tgh31 = Tg_r(3,1);
         Tgh32 = Tg_r(3,2);
         Tgh33 = Tg_r(3,3);
      elseif (i2 == 2*Neb+1)
         ibl = 3;
         Tgh21 = -swt(:,ibl)*Tg_r(1,1) + cwt(:,ibl)*Tg_r(2,1);
         Tgh22 = -swt(:,ibl)*Tg_r(1,2) + cwt(:,ibl)*Tg_r(2,2);
         Tgh23 = -swt(:,ibl)*Tg_r(1,3) + cwt(:,ibl)*Tg_r(2,3);
         Tgh31 = Tg_r(3,1);
         Tgh32 = Tg_r(3,2);
         Tgh33 = Tg_r(3,3);
      end

      s2_g = sgr_g                           ...
           + (Tg_r.')*[r(i2)*(cwt(:,ibl).'); ...
                      r(i2)*(swt(:,ibl).');  ...
                      zeros(1,N)]            ...
           + sxy_g;

      s12 = (s2_g - s1_g).';
      s = sqrt(s12(:,1).^2 + s12(:,2).^2 + s12(:,3).^2);
      sLu = s/oneLu;

      Qss = twosiggam*((0.5*sLu).^third).*besselk(third,sLu);
      dQss = -(twosiggam/oneLu)*((0.5*sLu).^third).*besselk(-2*third,sLu);
      dQsss = dQss./s;

      % Intermediate correlations in global coordinates.
      Qxx = Qss + 0.5*s.*dQss - 0.5*(s12(:,1).^2).*dQsss;
      Qxy = -0.5*s12(:,1).*s12(:,2).*dQsss;
      Qyx = Qxy;
      Qyy = Qss + 0.5*s.*dQss - 0.5*(s12(:,2).^2).*dQsss;
      Qzx = -0.5*s12(:,3).*s12(:,1).*dQsss;
      Qzy = -0.5*s12(:,3).*s12(:,2).*dQsss;

      % Transform the blade components from the global to the
      % uz,ut (hub z,y) frame.
      Q(:,icol+1) = Tgh31.*Qxx + Tgh32.*Qyx + Tgh33.*Qzx; % ux,uz
      Q(:,icol+2) = Tgh21.*Qxx + Tgh22.*Qyx + Tgh23.*Qzx; % ux,ut
      Q(:,icol+3) = Tgh31.*Qxy + Tgh32.*Qyy + Tgh33.*Qzy; % uy,uz
      Q(:,icol+4) = Tgh21.*Qxy + Tgh22.*Qyy + Tgh23.*Qzy; % uy,ut

   end

end

ival0 = 8*Nut*NN;

for i1 = 1:Nut

   ival = ival0 + 4*Nut*(i1-1);

   s1_g = [0;0;zut_g(i1)];

   for i2 = 1:Nut

      icol = ival + 4*(i2-1);

      s2_g = [0;0;zut_g(i2)] + sxy_g;
      s12 = (s2_g - s1_g).';

      if (i1 == i2)
         s12(1) = 1e-3*s12(2);  % Prevent divide-by-zero error.
      end

      s = sqrt(s12(:,1).^2 + s12(:,2).^2 + s12(:,3).^2);
      sLu = s/oneLu;

      Qss = twosiggam*((0.5*sLu).^third).*besselk(third,sLu);
      dQss = -(twosiggam/oneLu)*((0.5*sLu).^third).*besselk(-2*third,sLu);
      dQsss = dQss./s;

      Q(:,icol+1) = Qss + 0.5*s.*dQss - 0.5*(s12(:,1).^2).*dQsss;  % ux,ux
      Q(:,icol+2) = -0.5*s12(:,1).*s12(:,2).*dQsss;                % ux,uy
      Q(:,icol+3) = Q(:,icol+2);                                   % uy,ux
      Q(:,icol+4) = Qss + 0.5*s.*dQss - 0.5*(s12(:,2).^2).*dQsss;  % uy,uy

   end

end


