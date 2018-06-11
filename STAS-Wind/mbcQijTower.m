function Qpsi = mbcQijTower (Neb,Nut,Q,nf,df,psi0,Omega)
%
% The cross-correlations between blades and tower transform as
% Q^psi = T_B,p^psi(0) Q_ij(s_pq(tau),0)
% or
% Q^psi = Q_ij(s_pq(tau),0) [T_B,q(tau)]^T.
% That is, the transform appears only on the side of the
% correlation whose velocity components are representing the
% blades.
% 
% The input Q matrix is of the form
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
% Q[uy(elNut,tow),uy(elNut,tow),tau].
% The output matrix is identical, if we replace "bl1" with "0",
% "bl2" with "cos", and "bl3" with "sin". 
%
% The dimension is 4*nf-by-4*(2*Nut*3*Neb + Nut^2).
%
% Version:        Changes:
% --------        -------------
% 06.10.2015      Original code.
% 03.08.2017      Adapted for complex step derivatives.
%
% Version:        Verification:
% --------        -------------
% 06.10.2015      Memo AN 15.12.68
% 03.08.2017      
%
% Inputs:
% -------
% Neb             : Number of elements per blade.
% Q               : Correlation matrix.
% nf              : Number of frequencies for which an accurate
%                   result is desired.
% df              : Frequency bin width.
% psi0            : Initial rotor angle.
% Omega           : Rotational speed.
%
% Outputs:
% --------
% Qpsi            : Correlation matrix in MBC.

N = 4*nf;
NN = 3*Neb;
T0 = 1/df;
dt = T0/N;
tau = [[0:dt:T0/2].';[-(T0/2-dt):dt:-dt].']; % [0:dt:T0/2].'; %
Wt = Omega*tau;
twopi3 = 2*pi/3;
fourpi3 = 4*pi/3;
c0 = cos(psi0);
s0 = sin(psi0);
c2 = cos(twopi3 + psi0);
s2 = sin(twopi3 + psi0);
c4 = cos(fourpi3 + psi0);
s4 = sin(fourpi3 + psi0);
cp = zeros(N,3); % zeros(N/2+1,3);
sp = zeros(N,3); % zeros(N/2+1,3);
cwt = cos(Wt);
swt = sin(Wt);
cp(:,1) = cwt*c0 - swt*s0;
sp(:,1) = swt*c0 + cwt*s0;
cp(:,2) = cwt*c2 - swt*s2;
sp(:,2) = swt*c2 + cwt*s2;
cp(:,3) = cwt*c4 - swt*s4;
sp(:,3) = swt*c4 + cwt*s4;

Qpsi = zeros(N,4*(6*Nut*Neb + Nut^2));

for im1 = 1:3  % MBC 0, cos, sin.
   ival1 = 4*Neb*Nut*(im1-1);

   if (im1 == 1)
      v1 = [1 1 1];
   elseif (im1 == 2)
      % Get the cosines at tau = 0.
      v1 = 2*[cp(1,1) cp(1,2) cp(1,3)];
   else
      % Get the sines at tau = 0.
      v1 = 2*[sp(1,1) sp(1,2) sp(1,3)];
   end

   for i1 = 1:Neb
      ival2 = 4*Nut*(i1-1);
      
      for i2 = 1:Nut
         ival3 = 4*(i2-1);

         icol = [ival2+ival3           ...
                 ival2+ival3+4*Neb*Nut ...
                 ival2+ival3+8*Neb*Nut];    % The ref columns of Q.
         ipsi = ival1 + ival2 + ival3;      % The ref column of Qpsi.

         for j = 1:4
            Qpsi(:,ipsi+j) = (v1(1)*Q(:,icol(1)+j) ...
                            + v1(2)*Q(:,icol(2)+j) ...
                            + v1(3)*Q(:,icol(3)+j))/3;
         end

      end

   end

end

ival0 = 4*Nut*NN;

for i1 = 1:Nut
   ival1 = 4*NN*(i1-1);

   for im1 = 1:3
      ival2 = 4*Neb*(im1-1);

      if (im1 == 1)
         v2 = ones(N,3);
      elseif (im1 == 2)
         % Get the cosines.
         v2 = 2*cp;
      else
         % Get the sines.
         v2 = 2*sp;
      end

      for i2 = 1:Neb
         ival3 = 4*(i2-1);

         icol = [ival0+ival1+ival3       ...
                 ival0+ival1+ival3+4*Neb ...
                 ival0+ival1+ival3+8*Neb];
         ipsi = ival0 + ival1 + ival2 + ival3;

         for j = 1:4
            Qpsi(:,ipsi+j) = (Q(:,icol(1)+j).*v2(:,1) ...
                            + Q(:,icol(2)+j).*v2(:,2) ...
                            + Q(:,icol(3)+j).*v2(:,3))/3;
         end

      end

   end

end

ival0 = 8*Nut*NN;

Qpsi(:,ival0+1:ival0+4*(Nut^2)) = Q(:,ival0+1:ival0+4*(Nut^2));

