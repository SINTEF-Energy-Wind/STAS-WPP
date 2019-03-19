clear;

load 'LTMnorms.txt';
[length,time,mass,current,voltage,        ...
 velocity,force,power,stress,ndens,nvisc, ...
 stiffness,damping,resistance,inductance, ...
 capacitance,flux] = getLTMnorms ('LTMnorms.txt');

Neta = 84;
Naero = 126;

fmax = 2         * time;
vs0 = [33000 0].'  / voltage;

% Load closed-loop files.
nm = 'V117';
eval(["load 'Lpsi_" nm ".bin';"]);
eval(["load 'Apsi_" nm ".bin';"]);
eval(["load 'Bpsi_" nm ".bin';"]);
eval(["load 'Cpsi_" nm ".bin';"]);
eval(["load 'Dpsi_" nm ".bin';"]);
eval(["load 'dret_" nm ".bin';"]);
eval(["load 'xpsi_" nm ".txt';"]);
eval(["xpsi = xpsi_" nm ";"]);

cstate = [1 2 3 4 5 6 7 8 9 10 11 12 15 16 17 18].';
V0 = 12    / velocity;
Vref = 11.7 / velocity;
Pc0 = 1e7  / power;
Pcref = 1e7  / power;

% Partition out deleted DOFs.
Lr = Lpsi(dret,dret);
Ar = Apsi(dret,dret);
Br = Bpsi(dret,:);
Cr = Cpsi(:,dret);
DD = Dpsi;
xr = xpsi(dret);

%clear Lpsi Apsi Bpsi Cpsi Dpsi xpsi;

% Partition out control DOFs.
inds = [1:size(dret,1)].';
ind = inds(dret > 2*Neta+Naero+8+25);
Ncdof = size(ind,1);

[Lp,ret,cr] = partitionMatrix (Lr,ind,ind);
[Ap,rr,cr]  = partitionMatrix (Ar,ind,ind);
[Bp,rr,cr]  = partitionMatrix (Br,ind,[]);
[Cp,rr,cr]  = partitionMatrix (Cr,[],ind);
[xp,rr,cr]  = partitionMatrix (xr,ind,[]);
Nrdof = size(ret,1);

LL = Lp(1:Nrdof,1:Nrdof);
AA = Ap(1:Nrdof,1:Nrdof);
Bc = Ap(1:Nrdof,Nrdof+[1:Ncdof]);
B1 = Bp(1:Nrdof,:);
CC = Cp(:,1:Nrdof);
BB = [B1 Bc];
xs0 = xp(1:Nrdof);
xc0 = xp(Nrdof+[1:Ncdof]);

Nu1 = size(B1,2);

clear Lp Ap Bp Cp LAp LBp B1;

% =============================================================
% Do a modal reduction, and transform to a real state space.
LA = LL\AA;
LB = LL\BB;
[slap,shp,ifrq] = eigVal (LA);

% Some of the very closely spaced eigenfrequencies may lie out of order.
% This causes problems when I assume that conjugate modes are mirrored
% in the sorted mode shape matrix.  The solution is to force complex
% conjugacy by copying half the complex mode shapes and eigenvalues to
% the other half.  This is done in the 'for' loop below.
Phi = shp(:,ifrq);
Psi = inv(Phi);

iY = speye(Nrdof);
for im1 = 1:floor(Nrdof/2)

   imn = Nrdof - (im1-1);

   if (abs(imag(slap(im1))) > 0)
      Phi(:,imn) = conj(Phi(:,im1));
      Psi(imn,:) = conj(Psi(im1,:));
      slap(imn) = conj(slap(im1));
      iY([im1 imn],[im1 imn]) = [1 1;-i i];
   end

end

indfs = [1:Nrdof].';
idyn = indfs(abs(slap)<=fmax*2*pi);

ii = [1:Nrdof].';
jj = ii;
Lambda = sparse(ii,jj,slap);

Y = inv(iY);
Lamst = real(iY*Lambda*Y);  % will be real, but force real to eliminate noise.
YPB = real(iY*Psi*LB);

% =============================================================
% Set up the nonlinear controller.
c = STASControl_DTU10MW ();

% =============================================================
% Time simulation of the linked equations.

% Reduce to rotor-average wind speed and control state inputs.
iV = 618+[1:3:3*48-2];
ic = 618 + 3*48 + 6;
Bu = [sum(YPB(:,iV),2) YPB(:,ic+[1:Ncdof])];

inds = [1:size(dret,1)].';
iW  = inds(dret==164);
ib  = inds(dret==82);
iis = [inds(dret==294+8+17);inds(dret==294+8+18)].';

dt = 1/(fmax*50);
ts = [-30:dt*time:90].'   / time;
Nt = size(ts,1);

xc = xc0;

dzs = zeros(Nrdof,1);

% Extract the quasi-static part.
[jnk,iqs,cr] = partitionMatrix (indfs,idyn,[]);
Lamqs = Lamst(iqs,iqs);
LqsB = -Lamqs\Bu(iqs,:);

return

fid = fopen('out.txt','w');
for it = 1:Nt

   t = ts(it);

%   VV   = V0;
%   VV1  = V0;
   if (t < 30/time)
      VV   = V0;
      VV1  = V0;
   else
      VV   = V0 - 1/velocity;
      VV1  = V0 - 1/velocity;
   end

   PPc  = Pc0;
   PPc1 = Pc0;
%   if (t < 30/time)
%      PPc   = Pc0;
%      PPc1  = Pc0;
%   else
%      PPc   = Pc0 + 1.e6/power;
%      PPc1  = Pc0 + 1.e6/power;
%      PPc  = Pc0 + min( (ts(it-1)- 30)/5 , 1 )*(1e6/power);
%      PPc1 = Pc0 + min( (ts(it)  - 30)/5 , 1 )*(1e6/power);
%   end

   dV  = VV - Vref;
   dV1 = VV1 - Vref;
   dPc = PPc - Pcref;
   dPc1 = PPc1 - Pcref;   

   dxc = xc - xc0;
   dzs(iqs) = LqsB*[dV;dxc];
   k1s = Lamst(idyn,idyn)*dzs(idyn) + Bu(idyn,:)*[dV;dxc];
   xs = xs0 + real(Phi*Y)*dzs;
   W = xs(iW);
   b = xs(ib);
   isd = xs(iis(1));
   isq = xs(iis(2));
   Pe = isd*vs0(1) + isq*vs0(2);
   
   uc = [PPc;0;W;b;b;b;Pe;0;0;0;0;0;0];
   xcf = zeros(31,1);
   xcf(cstate) = xc;
   [k1c,ycout,jnkA,jnkB,jnkC,jnkD,jnk1,jnk2] =                        ...
         turbineControl (0,xcf,uc,c.cpar,c.cpct,                      ...
                         c.KeTab,c.WVTab,c.WPTab,c.bminTab,c.KTables, ...
                         c.KFTab,c.KSTab,c.KSqTab,c.KpiTab,c.KiiTab);
   VV2 = 0.5*(VV+VV1);
   PPc2 = 0.5*(PPc+PPc1);
   dV2 = VV2 - Vref;
   dPc2 = PPc2 - Pcref; 
   xc1 = xc + 0.5*k1c(cstate)*dt;
%xc1(2) = V0;
   dx1c = xc1 - xc0;
   dzs1 = dzs;
   dzs1(idyn) = dzs1(idyn) + 0.5*k1s*dt;
   dzs1(iqs)  = LqsB*[dV2;dx1c];
   k2s = Lamst(idyn,idyn)*dzs1(idyn) + Bu(idyn,:)*[dV2;dx1c];
   xs1 = xs0 + real(Phi*Y)*dzs1;
   W1 = xs1(iW);
   b1 = xs1(ib);
   isd1 = xs1(iis(1));
   isq1 = xs1(iis(2));
   Pe1 = isd1*vs0(1) + isq1*vs0(2);
   uc1 = [PPc2;0;W1;b1;b1;b1;Pe1;0;0;0;0;0;0];
   xcf1 = zeros(31,1);
   xcf1(cstate) = xc1;
   [k2c,ycout,jnkA,jnkB,jnkC,jnkD,jnk1,jnk2] =                        ...
         turbineControl (0,xcf1,uc1,c.cpar,c.cpct,                    ...
                         c.KeTab,c.WVTab,c.WPTab,c.bminTab,c.KTables, ...
                         c.KFTab,c.KSTab,c.KSqTab,c.KpiTab,c.KiiTab);
   xc = xc + k2c(cstate)*dt;
%xc(2) = V0;
   dzs(idyn) = dzs(idyn) + k2s*dt;



if (mod(it,10) == 1)
printf('%9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n', ...
       t,W,b,Pe,xc(2),xc(3),xc(7),xc(12));
fflush(stdout);
end
fprintf(fid,'%+5.6e %+5.6e %+5.6e %+5.6e %+5.6e %+5.6e\n', ...
        t,W,xc(3),b*180/pi,Pe,xc(2));


end

fclose('all');





