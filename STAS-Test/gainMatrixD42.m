
clear;

pkg load control;

Vmag = 10;
dt = 0.05;
fco = 0.5;

txt = '_P060_';
if (Vmag >= 10)
   Vstr = ['V' int2str(round(10*Vmag))];
else
   Vstr = ['V0' int2str(round(10*Vmag))];
end

eval(["load 'ADT" txt Vstr ".bin';"]);
eval(["load 'BDT" txt Vstr ".bin';"]);
eval(["load 'Ao" txt Vstr ".bin';"]);
eval(["load 'Bo" txt Vstr ".bin';"]);
eval(["load 'Co" txt Vstr ".bin';"]);
eval(["load 'Do" txt Vstr ".bin';"]);

Nx = size(ADT,1);
Ny = size(Co,1);
Nu = size(BDT,2);

% u         : 1: P, 2: beta, 3: yaw, 4: V, 5: thV, 6: Fw.
% y(sensors): 1: W, 2: beta0, 3: yaw, 4: Pe, 5-6: vnac, 7: V, 8: thV,
% y(other)  : 9: a, 10: FT, 11: Pa, 12-23: qfnd, 24-59: qbl, 60-61: isd,q.
iuc = [1:3].';
iue = [4:6].';
iyc = [4, 2, 3].';
iys = [1, 5, 6, 7, 8].';
Nys = size(iys,1);
Nyc = size(iyc,1);
Nuc = size(iuc,1);
Nue = size(iue,1);

Qdiag = [1.0, 0.1, 1.0];
Rdiag = [0.00005, 0.001, 0.001, 1.0, 0.0001];
QQ = diag(Qdiag);
RR = diag(Rdiag);
tol = 1.e-3;
G0 = zeros(Nx,Nys);
GG = RiccatiDT (2,ADT,BDT(:,iue),Co(iys,:),QQ,RR,G0,tol);

IGC = speye(Nx) - GG*Co(iys,:);
AG = IGC*ADT;
BG = BDT(:,iuc);  % The external commands.
CG = Co;
DG = Do(:,iuc);

% Check some step functions.
Nt = 1200;
yy = zeros(Nt,Nys*Nys);
yu = zeros(Nt,Nys*Nuc);
xy = zeros(Nx,Nys);
ys = zeros(Nys,1);
xu = zeros(Nx,Nuc);
u = zeros(Nuc,1);
ts = zeros(Nt,1);
for it = 2:Nt

   ts(it) = ts(it-1) + dt;

   if (it > 10)
      u = ones(Nuc,1);
      ys = ones(Nys,1);
   end
   for iu = 1:Nuc
      icn = Nys*(iu-1);
      u1 = zeros(Nuc,1);
      u1(iu) = 1;
      xu(:,iu) = AG*xu(:,iu) + BG*u1;
      yu(it,icn+[1:Nys]) = (CG(iys,:)*xu(:,iu) + DG(iys,:)*u1).';
   end
   for iy = 1:Nys
      icn = Nys*(iy-1);
      y1 = zeros(Nys,1);
      y1(iy) = ys(iy);
      xy(:,iy) = AG*xy(:,iy) + GG*y1;
      yy(it,icn+[1:Nys]) = (CG(iys,:)*xy(:,iy)).';
   end

end

eval(["save ('-binary','AG" txt Vstr ".bin','AG');"]);
eval(["save ('-binary','BG" txt Vstr ".bin','BG');"]);
eval(["save ('-binary','CG" txt Vstr ".bin','CG');"]);
eval(["save ('-binary','DG" txt Vstr ".bin','DG');"]);
eval(["save ('-binary','GG" txt Vstr ".bin','GG');"]);

% y(sensors): 1: W, 2: beta0, 3: yaw, 4: Pe, 5-6: vnac, 7: V, 8: thV,
clrz = [0.0 0.0 0.0; ...
        0.0 0.0 0.4; ...
        0.0 0.0 0.7; ...
        0.4 0.0 0.0; ...
        0.7 0.0 0.0];

figure(1);
clf;
for iyin = 1:Nys
   icn = Nys*(iyin-1);
   subplot(2,3,iyin);
   hold on;
   for iyout = 1:Nys
      plot (ts,yy(:,icn+iyout),'color',clrz(iyout,:),'linewidth',2.0);
   end
   hold off;
end

figure(2);
clf;
for iuin = 1:Nuc
   icn = Nys*(iuin-1);
   subplot(1,3,iuin);
   hold on;
   for iyout = 1:Nys
      plot (ts,yu(:,icn+iyout),'color',clrz(iyout,:),'linewidth',2.0);
   end
   hold off;
end

