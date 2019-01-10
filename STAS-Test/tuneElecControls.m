% Turbine electrical system feeding a static capacitor to ground, and series
% resistance and inductance connected to a voltage source.

clear;

epar = STASElectric_IEA10MW ();

ioff = [11 13 23 37 41];

np   = epar(1);
Ig   = epar(2:5);
Lg   = epar(6:9);
Rg   = epar(10);
lamr = epar(11);
Cdc  = epar(ioff(1)+1);
etac = epar(ioff(1)+2);
arat = epar(ioff(2)+1);
It   = epar(ioff(2)+[2:5]);
Lt   = epar(ioff(2)+[6:9]);
Rt   = epar(ioff(2)+10);
Ige  = epar(ioff(3)+[1:4]);
Lge  = epar(ioff(3)+[5:8]);
lamre= epar(ioff(3)+9);
KP   = epar(ioff(3)+10);
KI   = epar(ioff(3)+11);
KF   = epar(ioff(3)+12);
ag   = epar(ioff(3)+13);
aw   = epar(ioff(3)+14);
av   = epar(ioff(4)+1);
KPe  = epar(ioff(4)+2);
KIe  = epar(ioff(4)+3);
weh  = epar(ioff(4)+4);
ap   = epar(ioff(5)+1);
ais  = epar(ioff(5)+2);
avs  = epar(ioff(5)+3);
adc  = epar(ioff(5)+4);
KpDC = epar(ioff(5)+5);
KiDC = epar(ioff(5)+6);
KpQ  = epar(ioff(5)+7);
KiQ  = epar(ioff(5)+8);
Kppd = epar(ioff(5)+9);
Kipd = epar(ioff(5)+10);
Kppq = epar(ioff(5)+11);
Kipq = epar(ioff(5)+12);
KFp  = epar(ioff(5)+13);
Ite  = epar(ioff(5)+[14:17]);
Lte  = epar(ioff(5)+[18:21]);

Cn = 0;
Ln = 1*Lt(2)/(arat^2);
Rn = 0.5*Rt/(arat^2);

Nxs  = [2 1 2 5 4 11].';

wg   = 0.9089*0.5*np;
we   = weh;
th_e = 0;
vn   = [3933/arat;0];
ihg  = [1696;-4953];
Vhdc = 6500;
Qh   = 0;

spn = [0 -1;1 0];
I2 = eye(2);
wLR = we*Ln*spn + Rn*I2;
IwLRC = I2 + wLR*we*Cn*spn;

% Initialize the residual.
x0 = [ihg;Vhdc;ihg;ihg;0;0;we;th_e;vn;0;ihg;arat*ihg;vn;Vhdc;0;0;0;0];
%load 'x.txt';
%x0 = x;

vs = IwLRC\(vn + wLR*arat*x0(4:5));
th_e = atan2c(vs(2),vs(1));
yin = [wg;we;th_e;vs;ihg;Vhdc;Qh];
[dxdt,yout,A,B,C,D] = buildTurbineElectric (0,x0,yin,epar);

Nx  = size(A,1);
Ny  = size(C,1);
Nu  = size(B,2);

% Inactive states.
%deldofs = [sum(Nxs(1:4))+4 ...   % PLL Psie.
%           ].';
deldofs = [];
Nret = Nx - size(deldofs,1);
[jnk,ret,jnkc] = partitionMatrix ([1:Nx].',deldofs,[]);
x    = x0;
Res  = -dxdt(ret);
Rval = (Res.')*Res;

cnv   = eps^(0.2);
Ns    = 30;
bta   = [1.0*ones(1,100)].';

conv  = 0;
iter  = 0;
Niter = 0;
while ((Rval > cnv) && (iter < Ns))
   iter = iter + 1;

iter
Rval

   % Compute the tangent dynamics at the latest point.
   vs = IwLRC\(vn + wLR*arat*x(4:5));
   th_e = atan2c(vs(2),vs(1));
   yin = [wg;we;th_e;vs;ihg;Vhdc;Qh];
   [dxdt,yout,A,B,C,D] = buildTurbineElectric (1,x,yin,epar);

   % Modify the matrices to include the vs <=> (vn, is) relationship.
   bb = B(:,4:5)*(IwLRC\eye(2));
   A(:,4:5) = A(:,4:5) + B(:,4:5)*(IwLRC\(wLR*arat));
   B(:,4:5) = bb;  % Replace vs with vn as the global input.

   dRdxn = -A(ret,ret);

   % ----------------------------------------------------------------
   % Apply Newton's method on the MBC-transformed equations.
   lam = 1;
   dx = -dRdxn\Res;
   lflg = 0;
   litmax = 20;
   liter = 0;
   while ((lflg == 0) && (liter < litmax))
      liter = liter + 1;
%'-----'
%liter

      xn      = x0;
      xn(ret) = x(ret) + bta(iter)*lam*dx;

      vs = IwLRC\(vn + wLR*arat*xn(4:5));
      th_e = atan2c(vs(2),vs(1));
      yin = [wg;we;th_e;vs;ihg;Vhdc;Qh];
      [dxndt,ynout,jnkA,jnkB,jnkC,jnkD] = buildTurbineElectric (0,xn,yin,epar);
      Res   = -dxndt(ret);
      R1    = (Res.')*Res;

%[dx x(ret) xn(ret)]
%[dxdt(ret) dxndt(ret) Res]
%[R1 Rval]

      if (R1 < Rval)
         % OK!  Prepare for the next iteration.
         lflg  = 1;
         Rval  = R1;
         x     = xn;
         dxdt  = dxndt;
      else
         % Backtrack.
         lam = 0.5*lam;
         if (liter == litmax)
            [iter R1 Rval]
            'Warning, proceeding without lambda convergence.'
            lflg = 1;
            x = xn;
            dxdt = dxndt;
return
         end
      end

   end  % Newton inner.

   if (iter == Ns) && (Rval > cnv)
      Rval
      'Warning, max iterations, proceeding without Rval convergence.'
return
   end

end  % Newton outer.

Rval

save('-ascii','x.txt','x');
vs = IwLRC\(vn + wLR*arat*x(4:5));
th_e = atan2c(vs(2),vs(1));
yin = [wg;we;th_e;vs;ihg;Vhdc;Qh];
[dxdt,yout,A,B,C,D] = buildTurbineElectric (1,x,yin,epar);

bb = B(:,4:5)*(IwLRC\eye(2));
A(:,4:5) = A(:,4:5) + B(:,4:5)*(IwLRC\(wLR*arat));
B(:,4:5) = bb;  % Replace vs with vn as the global input.

A = A(ret,ret);
B = B(ret,:);
C = C(:,ret);

[slap,shp,ifrq] = eigVal(A);



freqs = [[0.001:0.001:0.1] [0.11:0.01:1] [1.1:0.1:100] [101:1:1000]].';
Nf = size(freqs,1);
dVdcdVhdc = zeros(Nf,1);
dQdQh = zeros(Nf,1);
%digdihg = zeros(Nf,2);
%dipdihp = zeros(Nf,2);
%dwemdwe = zeros(Nf,1);

for ifreq = 1:Nf

   freq = freqs(ifreq);
   omega = 2*pi*freq;

   dxdu = (i*omega*eye(Nret) - A)\B;

   dVdcdVhdc(ifreq) = dxdu(3,8);
%   digdihg(ifreq,:) = [dxdu(1,6) dxdu(2,7)];

   disdQh = arat*dxdu(4:5,9);
   dvsdQh = IwLRC\(wLR*arat*dxdu(4:5,9));
   dQdQh(ifreq) = dvsdQh(2)*arat*x(4) - dvsdQh(1)*arat*x(5) + vs(2)*disdQh(1) - vs(1)*disdQh(2);

end


figure(1);
hold on;
semilogx(freqs,sqrt(real(dVdcdVhdc).^2 + imag(dVdcdVhdc).^2));
semilogx(freqs,sqrt(real(dQdQh).^2 + imag(dQdQh).^2));
figure(2);
hold on;
semilogx(freqs,atan2(imag(dVdcdVhdc),real(dVdcdVhdc))/pi);
semilogx(freqs,atan2(imag(dQdQh),real(dQdQh))/pi);
