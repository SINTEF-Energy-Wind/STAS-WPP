function Qpsi = mbcQij (Neb,Q,nf,df,psi0,Omega)
%
% The correlation matrix of turbulence transforms as
% Q^psi = T_B,p^psi(0) Q_ij(s_pq(tau),0) [T_B,q(tau)]^T.
% i and j denote the z,t components, p and q denote the points
% between which the correlation is being taken.
%
% The input Q matrix is of the form
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
% Q[ut(elNeb,bl1),ut(elNeb,bl3),tau].
% This is a compressed format which takes advantage of the
% assumed isotropic turbulence such that the correlation
% between blade 2 and blade 3 is the same as that between
% blade 1 and blade 2, and so on.
%
% The output Qpsi matrix is of the form
% Q[uz(el1  ,0  ),uz(el1  ,0  ),tau]
% Q[uz(el1  ,0  ),ut(el1  ,0  ),tau]
% Q[ut(el1  ,0  ),uz(el1  ,0  ),tau]
% Q[ut(el1  ,0  ),ut(el1  ,0  ),tau]
% Q[uz(el1  ,0  ),uz(el2  ,0  ),tau]
%               ...
% Q[ut(el1  ,0  ),ut(elNeb,0  ),tau]
% Q[uz(el1  ,0  ),uz(el1  ,cos),tau]
%               ...
% Q[ut(el1  ,0  ),ut(elNeb,sin),tau]
% Q[uz(el2  ,0  ),uz(el1  ,0  ),tau]
%               ...
% Q[ut(elNeb,0  ),ut(elNeb,sin),tau]
% Q[uz(el1  ,cos),uz(el1  ,0  ),tau]
%               ...
% Q[ut(elNeb,sin),ut(elNeb,sin),tau].
% That is, it includes all the components: collective, cos, and
% sin.  But, these are each of length Neb, so the number of
% columns in Qpsi is 4*(3*Neb)^2.
%
% Version:        Changes:
% --------        -------------
% 15.02.2015      Original code.
% 04.08.2017      Adapted for complex step derivatives.
%
% Version:        Verification:
% --------        -------------
% 15.02.2015      Memo 15.12.19
% 04.08.2017      
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
T0 = 1/df;
dt = T0/N;
tau = [[0:dt:T0/2].';[-(T0/2-dt):dt:-dt].'];
Wt = Omega*tau;
twopi3 = 2*pi/3;
fourpi3 = 4*pi/3;
c0 = cos(psi0);
s0 = sin(psi0);
c2 = cos(twopi3 + psi0);
s2 = sin(twopi3 + psi0);
c4 = cos(fourpi3 + psi0);
s4 = sin(fourpi3 + psi0);
cp = zeros(N,3);
sp = zeros(N,3);
cwt = cos(Wt);
swt = sin(Wt);
cp(:,1) = cwt*c0 - swt*s0;
sp(:,1) = swt*c0 + cwt*s0;
cp(:,2) = cwt*c2 - swt*s2;
sp(:,2) = swt*c2 + cwt*s2;
cp(:,3) = cwt*c4 - swt*s4;
sp(:,3) = swt*c4 + cwt*s4;

Qpsi = zeros(N,36*(Neb^2));

% Q has been stored in a form that omits repeated components
% of the correlations.  The mapping from the values in the
% Q array to the full list of correlations between the three
% blades is
% Full:   11  12  13  21  22  23  31  32  33
% Stored: 11  12  13  13  11  12  12  13  11
for im1 = 1:3
   ival1 = 12*(Neb^2)*(im1-1);

   if (im1 == 1)
      v1 = [1 1 1];
   elseif (im1 == 2)
      % Get the cosines at tau = 0.
      v1 = 2*[cp(1,1) cp(1,2) cp(1,3)];
   else
      % Get the sines at tau = 0.
      v1 = 2*[sp(1,1) sp(1,2) sp(1,3)];
   end

   for ie1 = 1:Neb
      ival2 = 12*Neb*(ie1-1);
      for im2 = 1:3
         ival3 = 4*Neb*(im2-1);

         if (im2 == 1)
            v2 = ones(N,3);
         elseif (im2 == 2)
            % Get the cosines.
            v2 = 2*cp;
         else
            % Get the sines.
            v2 = 2*sp;
         end

         for ie2 = 1:Neb
            ival4 = 4*(ie2-1);
            ivb = [ival2+ival4 ival2+ival4+4*Neb ival2+ival4+8*Neb];
            ipsi = ival1 + ival2 + ival3 + ival4;
            icol = [ivb(1) ivb(2) ivb(3) ...
                    ivb(3) ivb(1) ivb(2) ...
                    ivb(2) ivb(3) ivb(1)];

            for j = 1:4

               Qpsi(:,ipsi+j) =                          ...
                          (v1(1)*Q(:,icol(1)+j).*v2(:,1) ...
                        +  v1(1)*Q(:,icol(2)+j).*v2(:,2) ...
                        +  v1(1)*Q(:,icol(3)+j).*v2(:,3) ...
                        +  v1(2)*Q(:,icol(4)+j).*v2(:,1) ...
                        +  v1(2)*Q(:,icol(5)+j).*v2(:,2) ...
                        +  v1(2)*Q(:,icol(6)+j).*v2(:,3) ...
                        +  v1(3)*Q(:,icol(7)+j).*v2(:,1) ...
                        +  v1(3)*Q(:,icol(8)+j).*v2(:,2) ...
                        +  v1(3)*Q(:,icol(9)+j).*v2(:,3))/9;

            end
         end
      end
   end
end



