function [dxdt,A,By] = Vidyn (psiFlag,x,Vzts,Viy,r,Do,azi,Omega)
%
% State equations for the dynamic induced velocity, with the option
% to use multi-blade coordinates.
%
%   States:           y vector:         u vector:
%   Vih   1:2         Vzts    1:3
%   Vi    3:4         Viy     4:5
%                     r       6
% (MBC: repeat the above for each blade.)
%                     azi     7 or 19
%                     Omega   8 or 20
%                     D       9 or 21-23 
%
% Version:        Changes:
% --------        -------------
% 02.01.2018      Original code.
%
% Version:        Verification:
% --------        -------------
% 02.01.2018      Checked A, By by complex step using dx/dt output.
%                 dx/dt matches BEMNL.m, checked for both psiFlag=0,1.
%
% Inputs:
% -------
% psiFlag         : 0: independent elements.  1: MBC coordinates.
% x               : States Vih_z,t and Vi_z,t.
% Vzts            : Local flow in rotorplane z,t,s coordinates.
% Viy             : QS induced velocity including yaw correction.
% xer             : Projected x,y,z coordinates of the element on the
%                   rotorplane.
% r,Do            : Representative values for computing the time
%                   constants.
% azi             : Rotor azimuth.
% Omega           : Rotor speed.
%
% Outputs:
% --------
% dxdt            : LHS of state equations.
% A,By            : State space matrices.

if (psiFlag == 0)

   % Independent elements.  Assume input is for a single element.
   dxdt = zeros (4,1);
   A    = zeros (4,4);
   By   = zeros (4,9);

   iVzts  = 0;
   iViy   = 3;
   ir     = 5;
%   iazi   = 6;  % Unused.
%   iOmega = 7;  % Unused.
   iD     = 8;

   Vih = x(1:2);
   Vi  = x(3:4);

   Vmag = sqrt(Vzts(1)^2 + Vzts(2)^2 + Vzts(3)^2);
   aa   = -Vi(1)/Vzts(1);
   D2V  = Do/(2*Vmag);
   tau1 = (1.1/(1 - 0.3*aa))*D2V;
   tau2 = (0.39 - 0.26*((2*r/Do)^2))*tau1;
   it1 = 1/tau1;
   it2 = 1/tau2;

   dt1dV   = -tau1/Vmag;
   dt1dvz  =  (0.33/((1 - 0.3*aa)^2))*D2V*Vi(1)/(Vzts(1)^2);
   dt1dviz = -(0.33/((1 - 0.3*aa)^2))*D2V/Vzts(1);
   dt1dD   =  tau1/Do;
   dt2dt1  =  0.39 - 0.26*((2*r/Do)^2);
   dt2dr   = -0.26*8*r*tau1/(Do^2);
   dt2dD   =  0.26*8*(r^2)*tau1/(Do^3);
   dVdvz   =  Vzts(1)/Vmag;
   dVdvt   =  Vzts(2)/Vmag;
   dVdvs   =  Vzts(3)/Vmag;

   dxdt(1) = -it1*Vih(1) + 0.4*it1*Viy(1);
   dxdt(2) = -it1*Vih(2) + 0.4*it1*Viy(2);
   dxdt(3) =  it2*Vih(1) - it2*Vi(1) + 0.6*it2*Viy(1);
   dxdt(4) =  it2*Vih(2) - it2*Vi(2) + 0.6*it2*Viy(2);

   Amat  = [-it1 0;it2 -it2];
   At1   = [it1^2 0;0 0];
   dAdt2 = [0 0;-it2^2 it2^2];

   Bvec  = [0.4*it1;0.6*it2];
   Bt1   = [-0.4*it1^2;0];
   dBdt2 = [0;-0.6*it2^2];

   dAdt1  = At1 + dAdt2*dt2dt1;
   dBdt1  = Bt1 + dBdt2*dt2dt1;
   dAdV   = dAdt1*dt1dV;
   dBdV   = dBdt1*dt1dV;
   dAdvz  = dAdt1*dt1dvz;
   dBdvz  = dBdt1*dt1dvz;
   dAdviz = dAdt1*dt1dviz;
   dBdviz = dBdt1*dt1dviz;
   dAdr   = dAdt2*dt2dr;
   dAdD   = dAdt1*dt1dD + dAdt2*dt2dD;

   for icomp = 1:2

      irow = [icomp icomp+2].';
      A(irow,irow)        = Amat;
      By(irow,iViy+icomp) = Bvec;

      icol = 3;
      ind  = icomp;
      A(irow,icol) = A(irow,icol)                     ...
                   + dAdt1*dt1dviz*[Vih(ind);Vi(ind)] ...
                   + dBdt1*dt1dviz*Viy(ind);
         
      icol = iVzts + 1;
      By(irow,icol) = By(irow,icol)                   ...
                    + dAdt1*dt1dvz*[Vih(ind);Vi(ind)] ...
                    + dAdV *dVdvz *[Vih(ind);Vi(ind)] ...
                    + dBdt1*dt1dvz*Viy(ind)           ...
                    + dBdV *dVdvz *Viy(ind);

      icol = iVzts + 2;
      By(irow,icol) = By(irow,icol)                 ...
                    + dAdV*dVdvt*[Vih(ind);Vi(ind)] ...
                    + dBdV*dVdvt*Viy(ind);

      icol = iVzts + 3;
      By(irow,icol) = By(irow,icol)                 ...
                    + dAdV*dVdvs*[Vih(ind);Vi(ind)] ...
                    + dBdV*dVdvs*Viy(ind);

      icol = ir + 1;
      By(irow,icol) = By(irow,icol)                  ...
                    + dAdt2*dt2dr*[Vih(ind);Vi(ind)] ...
                    + dBdt2*dt2dr*Viy(ind);

      icol = iD + 1;
      By(irow,icol) = By(irow,icol)                  ...
                    + dAdt1*dt1dD*[Vih(ind);Vi(ind)] ...
                    + dAdt2*dt2dD*[Vih(ind);Vi(ind)] ...
                    + dBdt1*dt1dD*Viy(ind)           ...
                    + dBdt2*dt2dD*Viy(ind);

   end

elseif (psiFlag == 1)

   % MBC triplet of elements.
   dxdt = zeros (12,1);
   A    = zeros (12,12);
   By   = zeros (12,23);

   i2a = [1:2:5].';
   i2b = [2:2:6].';
   i3a = [1:3:7].';
   i3b = [2:3:8].';
   i3c = [3:3:9].';

   indx = [0 4 8].';
   indb = [0 2 4].';

   Vih = x([1 2 5 6 9 10]);
   Vi  = x([3 4 7 8 11 12]);

   Vmag = sqrt(Vzts(i3a).^2 + Vzts(i3b).^2 + Vzts(i3c).^2);
   aa   = -Vi(i2a)./Vzts(i3a);
   D2V  = Do./(2*Vmag);

%'NL dVdz'
%[real(Vmag) imag(Vmag)/sqrt(eps)]

   tau1 = (1.1./(1 - 0.3*aa)).*(Do./(2*Vmag));
   tau2 = (0.39 - 0.26*(2*r./Do).^2).*tau1;
   tm1 = mean(tau1)*ones(3,1);
   tm2 = mean(tau2)*ones(3,1);
   it1  = 1./tm1;
   it2  = 1./tm2;

%'NL dt1dz'
%[real(tau1) imag(tau1)/sqrt(eps)]

%'NL dtm1dz'
%[real(tm1) imag(tm1)/sqrt(eps)]

   caz = [cos(azi) cos(azi+(2*pi/3)) cos(azi+(4*pi/3))].';
   saz = [sin(azi) sin(azi+(2*pi/3)) sin(azi+(4*pi/3))].';

   psiB  = [1  caz(1)  saz(1); ...
            1  caz(2)  saz(2); ...
            1  caz(3)  saz(3)];

   dpsiB = [0 -saz(1)  caz(1); ...
            0 -saz(2)  caz(2); ...
            0 -saz(3)  caz(3)];

   Bpsi   = (1/3)*[    1         1         1;     ...
                    2*caz(1)  2*caz(2)  2*caz(3); ...
                    2*saz(1)  2*saz(2)  2*saz(3)];

   dBpsi  = (1/3)*[    0         0         0;     ...
                   -2*saz(1) -2*saz(2) -2*saz(3); ...
                    2*caz(1)  2*caz(2)  2*caz(3)];

   d2Bpsi = (1/3)*[    0         0         0;     ...
                   -2*caz(1) -2*caz(2) -2*caz(3); ...
                   -2*saz(1) -2*saz(2) -2*saz(3)];

   Atau = zeros (12,12);
   Atau(indx+1,indx+1) = diag(-it1);
   Atau(indx+3,indx+1) = diag(it2);
   Atau(indx+3,indx+3) = diag(-it2);
   Atau(indx+2,indx+2) = diag(-it1);
   Atau(indx+4,indx+2) = diag(it2);
   Atau(indx+4,indx+4) = diag(-it2);

   Btau = zeros (12,6);
   Btau(indx+1,indb+1) = diag(0.4*it1);
   Btau(indx+3,indb+1) = diag(0.6*it2);
   Btau(indx+2,indb+2) = diag(0.4*it1);
   Btau(indx+4,indb+2) = diag(0.6*it2);

   TpsiBx = zeros (12,12);
   TpsiBx(indx+1,indx+1) = psiB;
   TpsiBx(indx+2,indx+2) = psiB;
   TpsiBx(indx+3,indx+3) = psiB;
   TpsiBx(indx+4,indx+4) = psiB;

   TBpsix = zeros (12,12);
   TBpsix(indx+1,indx+1) = Bpsi;
   TBpsix(indx+2,indx+2) = Bpsi;
   TBpsix(indx+3,indx+3) = Bpsi;
   TBpsix(indx+4,indx+4) = Bpsi;

   dTpsiBx = zeros (12,12);
   dTpsiBx(indx+1,indx+1) = dpsiB;
   dTpsiBx(indx+2,indx+2) = dpsiB;
   dTpsiBx(indx+3,indx+3) = dpsiB;
   dTpsiBx(indx+4,indx+4) = dpsiB;

   dTBpsix = zeros (12,12);
   dTBpsix(indx+1,indx+1) = dBpsi;
   dTBpsix(indx+2,indx+2) = dBpsi;
   dTBpsix(indx+3,indx+3) = dBpsi;
   dTBpsix(indx+4,indx+4) = dBpsi;

   d2TBpsix = zeros (12,12);
   d2TBpsix(indx+1,indx+1) = d2Bpsi;
   d2TBpsix(indx+2,indx+2) = d2Bpsi;
   d2TBpsix(indx+3,indx+3) = d2Bpsi;
   d2TBpsix(indx+4,indx+4) = d2Bpsi;

   TBpsib = zeros (6,6);
   TBpsib(indb+1,indb+1) = Bpsi;
   TBpsib(indb+2,indb+2) = Bpsi;

   dTBpsib = zeros (6,6);
   dTBpsib(indb+1,indb+1) = dBpsi;
   dTBpsib(indb+2,indb+2) = dBpsi;

   Apsi = TpsiBx*(Atau*TBpsix - Omega*dTBpsix);
   Bmpsi = TpsiBx*Btau*TBpsib;

   dxdt = Apsi*x + Bmpsi*Viy;

   dtmdt = 1/3;

   dt1dV   = -(1.1./(1 - 0.3*aa)).*D2V./Vmag;
   dt1dvz  =  (0.33./((1 - 0.3*aa).^2)).*D2V.*Vi(i2a)./(Vzts(i3a).^2);
   dt1dviz = -(0.33./((1 - 0.3*aa).^2)).*D2V./Vzts(i3a);
   dt1dD   =  (1.1./(1 - 0.3*aa))./(2*Vmag);

   dt2dt1  =  (0.39 - 0.26*((2*r./Do).^2));
   dt2dr   = -0.26*8*r.*tau1./(Do.^2);
   dt2dD   =  0.26*8*(r.^2).*tau1./(Do.^3);

   dVdvz   =  Vzts(i3a)./Vmag;
   dVdvt   =  Vzts(i3b)./Vmag;
   dVdvs   =  Vzts(i3c)./Vmag;

%real(dVdvz)
%real(dVdvt)
%real(dVdvs)

   A(:,:) = Apsi;
   By(:,[4 5 10 11 16 17]) = Bmpsi;

   % dA_tau/dtau_mean and dB_tau/dtau_mean.  These may be computed
   % outside the loop, due to the averaging of the time constants.
   dAdtm1 = zeros (12,12);  % Direct derivative dAtau/dtm1
   dAdtm1(indx+1,indx+1) = diag(it1.^2);
   dAdtm1(indx+2,indx+2) = diag(it1.^2);

   dAdtm2 = zeros(12,12);
   dAdtm2(indx+3,indx+1) = diag(-it2.^2);
   dAdtm2(indx+3,indx+3) = diag(it2.^2);
   dAdtm2(indx+4,indx+2) = diag(-it2.^2);
   dAdtm2(indx+4,indx+4) = diag(it2.^2);

   dBdtm1 = zeros (12,6);
   dBdtm1(indx+1,indb+1) = diag(-0.4*it1.^2);
   dBdtm1(indx+2,indb+2) = diag(-0.4*it1.^2);

   dBdtm2 = zeros (12,6);
   dBdtm2(indx+3,indb+1) = diag(-0.6*it2.^2);
   dBdtm2(indx+4,indb+2) = diag(-0.6*it2.^2);

   At1   = dAdtm1*dtmdt;  % Same for all the tau's.
   dAdt2 = dAdtm2*dtmdt;
   Bt1   = dBdtm1*dtmdt;
   dBdt2 = dBdtm2*dtmdt;

   TAT2 = TpsiBx*dAdt2*TBpsix;
   TBT2 = TpsiBx*dBdt2*TBpsib;

   for it = 1:3

      ic6 = 6*(it-1);

      dAdt1 = At1 + dAdt2*dt2dt1(it);
      dBdt1 = Bt1 + dBdt2*dt2dt1(it);

      TAT1 = TpsiBx*dAdt1*TBpsix;
      TBT1 = TpsiBx*dBdt1*TBpsib;

      term_viz = TAT1*dt1dviz(it)*x + TBT1*dt1dviz(it)*Viy;
      term_vz  = TAT1*dt1dvz(it)*x          ...
               + TAT1*dt1dV(it)*dVdvz(it)*x ...
               + TBT1*dt1dvz(it)*Viy        ...
               + TBT1*dt1dV(it)*dVdvz(it)*Viy;
      term_vt  = TAT1*dt1dV(it)*dVdvt(it)*x ...
               + TBT1*dt1dV(it)*dVdvt(it)*Viy;
      term_vs  = TAT1*dt1dV(it)*dVdvs(it)*x ...
               + TBT1*dt1dV(it)*dVdvs(it)*Viy;
      term_r   = TAT2*dt2dr(it)*x + TBT2*dt2dr(it)*Viy;
      term_D   = TAT1*dt1dD(it)*x + TBT1*dt1dD(it)*Viy ...
               + TAT2*dt2dD(it)*x + TBT2*dt2dD(it)*Viy;

      icol = indx(it) + 3;
      A(:,icol)  = A(:,icol)  + term_viz;

      icol = ic6 + 1;
      By(:,icol) = By(:,icol) + term_vz; 

      icol = ic6 + 2;
      By(:,icol) = By(:,icol) + term_vt;

      icol = ic6 + 3;
      By(:,icol) = By(:,icol) + term_vs; 

      icol = ic6 + 6;
      By(:,icol) = By(:,icol) + term_r;

      icol = 20 + it;
      By(:,icol) = By(:,icol) + term_D;

   end

   icol = 19;  % azi.
   By(:,icol) = By(:,icol)                               ...
              - TpsiBx*(-dTBpsix*dxdt - Omega*d2TBpsix*x ...
              +         Atau*dTBpsix*x + Btau*dTBpsib*Viy);

   icol = 20;  % Omega.
   By(:,icol) = By(:,icol) - TpsiBx*dTBpsix*x;

end

