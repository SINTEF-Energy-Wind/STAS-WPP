% Run a set of steady-state calculations of generator efficiency, which is
% used in the windspeed observer.

clear;

[length,time,mass,current,voltage,        ...
 velocity,force,power,stress,ndens,nvisc, ...
 stiffness,damping,resistance,inductance, ...
 capacitance,flux] = getLTMnorms ('LTMnorms.txt');

nm = 'DTU10MW';

eval(['epar  = STASElectric_' nm ' ();']);
eval(['c     = STASControl_'  nm ' ();']);

ioff = [11 13 23 37 41 62];

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
Pnll = epar(ioff(6)+1);

Cpnom = 0.486;
Area  = pi*(c.Ro^2);
dens  = 1.225                     / ndens;
ihgqs = -[100:100:1000].'       / current;
wgs   = 0.5*np*0.628; % 0.5*np*[0.6:0.1:1].'      * time;

Ni = size(ihgqs,1);
Nw = size(wgs,1);

Pes  = zeros(Ni*Nw);
etas = zeros(Ni*Nw);

yin = zeros(9,1);
yin(2) = 50*2*pi                  * time;
yin(3) = 0;
yin(4) = 33000                    / voltage;
yin(5) = 0.1                      / voltage;
yin(6) = 0.1                      / current;
yin(8) = c.Vhdc;
yin(9) = 0                        / power;

x0 = zeros(25,1);
x0(1)  = yin(6);
x0(3)  = yin(8);
x0(4)  = yin(4)*arat;
x0(5)  = yin(5)*arat;
x0(6)  = yin(6);
x0(8)  = x0(6);
x0(11) = yin(3);
x0(12) = yin(4);
x0(13) = yin(5);
x0(14) = yin(4);
x0(15) = x0(4);
x0(16) = x0(5);
x0(17) = yin(4);
x0(18) = yin(5);
x0(19) = yin(8);
x0(20) = yin(5);
x0(21) = yin(8);
x0(22) = yin(8);
x0(23) = yin(9);
x0(24) = x0(15);
x0(25) = x0(16);

for iw = 1:Nw

   wg = wgs(iw);
   yin(1) = wg;
   x0(10) = wg;

   for ii = 1:Ni

      ind = Ni*(iw-1) + ii;

      ihgq = ihgqs(ii);
      yin(7) = ihgq;
      x0(2)  = yin(7);
      x0(7)  = yin(7);
      x0(9)  = yin(7);

      cnv = eps^0.6;
      Ns = 50;
      bta = ones(Ns,1);
      litmax = 20;
      [xs,dxs] = solveNewt (@(x) BTEfun(x,yin,epar),0,x0,cnv,Ns,bta,litmax);

      [dxdt,yout,A,B,C,D] = buildTurbineElectric (0,xs,yin,epar);
      Ps(ind) = (yout(2:3).')*yin(4:5) - Pnll;
      W = wg/(0.5*np);
      Tg = yout(1);
      etag(ind) = Ps(ind)/(Tg*W);

      printf('%10.4f %10.4f %10.4f %10.4f\n',W,Ps(ind),Tg*W,etag(ind));
fflush(stdout);

   end

end


