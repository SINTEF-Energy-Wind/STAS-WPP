function dxdt = RotCantTimestep (x,t,mes,kes,PB,Pns,WW0,adamp)
%
% States: q, dq/dt (each 6*Nel), azi.
%

Nx = size(x,1);
Ndof = (Nx - 1)/2;
Nel = Ndof/6;

% Some convenient indices.
i6a = [1:6:6*Nel-5].';
i6b = [2:6:6*Nel-4].';
i6c = [3:6:6*Nel-3].';
i6d = [4:6:6*Nel-2].';
i6e = [5:6:6*Nel-1].';
i6f = [6:6:6*Nel].';

% ===================================================================
% Define functions of time here...
F = zeros(Ndof,1);
if ((t > 0.099) && (t < 0.100))
   F(Ndof-6+2) = 1000*sin((t-0.099)*pi/0.001);
end
WW = min(WW0,(t/0.05)*WW0);  % Linear ramp up, > period of first mode.
% ===================================================================

r = sqrt(PB(1)^2 + PB(2)^2);

% Update qB and dqBdt.
azi = x(Nx);
ca = cos(azi);
sa = sin(azi);
TBg = [ca -sa 0; sa ca 0; 0 0 1];
qB = zeros(6,1);
qB(4:6) = thetaFromT (TBg);
qB(1:3) = TBg*[r;0;0];

dTBg = [-sa -ca 0; ca -sa 0; 0 0 0];
Gz = TBg*(dTBg.');
Gx = zeros(3,3);
Gy = zeros(3,3);

[jnk,dT] = dTdth (qB(4:6));
Fx = dT(:,1:3)*(TBg.');
Fy = dT(:,4:6)*(TBg.');
Fz = dT(:,7:9)*(TBg.');

GG = [-Gx(2,3) -Gy(2,3) -Gz(2,3); ...
       Gx(1,3)  Gy(1,3)  Gz(1,3); ...
      -Gx(1,2) -Gy(1,2) -Gz(1,2)];
FF = [-Fx(2,3) -Fy(2,3) -Fz(2,3); ...
       Fx(1,3)  Fy(1,3)  Fz(1,3); ...
      -Fx(1,2) -Fy(1,2) -Fz(1,2)];
dqBdt = zeros(6,1);
dqBdt(4:6) = -FF\GG*[0;0;WW];
dqBdt(1:3) = dTBg*[r;0;0]*WW;  % (Returns dqBdt(6) = WW, as expected.)

M = zeros(Ndof,Ndof);
G = zeros(Ndof,Ndof);
H = zeros(Ndof,1);
K = zeros(Ndof,1);
C = zeros(Ndof,Ndof);
for iel = 1:Nel

   ic6 = 6*(iel-1);

   if (iel == 1)
      rdof     = 0;
      qn1      = zeros(6,1);
      dqn1dt   = zeros(6,1);
      Pn1      = zeros(6,1);
      Pn1(4:6) = Pns(10:12);
   else
      rdof   = ic6 - 6;
      qn1    = x(ic6-6+[1:6]);
      dqn1dt = x(6*Nel+ic6-6+[1:6]);
      Pn1    = Pns(ic6-6+[1:6]);
   end

   qn2    = x(ic6+[1:6]);
   dqn2dt = x(6*Nel+ic6+[1:6]);
   Pn2    = Pns(ic6+[1:6]);

   dqdt = [dqBdt;dqn1dt;dqn2dt];

   Qu1  = Qunod (qn1,qB,Pn1,PB);
   Qu2  = Qunod (qn2,qB,Pn2,PB);
   dQu1 = dQudq (qn1,qB,Pn1,PB);
   dQu2 = dQudq (qn2,qB,Pn2,PB);

   [xe,TsB] = elementCSFromNodes (qn1,qn2,Pn1,Pn2);
   dTsB = derivElementCS (qn1,qn2,Pn1,Pn2,TsB);

   mu = getMu (qn1,qn2,Pn1,Pn2,TsB);
   dmu = dmudq (qn1,qn2,Pn1,Pn2,TsB,dTsB);

   [mme,gge,hhe,kke] = buildElementNL (mes,kes,dqdt,      ...
                                       Qu1,Qu2,dQu1,dQu2, ...
                                       mu,dmu,TsB,dTsB);

   cce = dampingC (adamp*kes,dmu);

   if (iel == 1)
      M(rdof+[1:6],rdof+[1:6]) = M(rdof+[1:6],rdof+[1:6]) + mme(13:18,13:18);
      G(rdof+[1:6],rdof+[1:6]) = G(rdof+[1:6],rdof+[1:6]) + gge(13:18,13:18);
      H(rdof+[1:6])            = H(rdof+[1:6])            + hhe(13:18);
      K(rdof+[1:6])            = K(rdof+[1:6])            + kke(13:18);
      C(rdof+[1:6],rdof+[1:6]) = C(rdof+[1:6],rdof+[1:6]) + cce(13:18,13:18);
   else
      M(rdof+[1:12],rdof+[1:12]) = M(rdof+[1:12],rdof+[1:12]) + mme(7:18,7:18);
      G(rdof+[1:12],rdof+[1:12]) = G(rdof+[1:12],rdof+[1:12]) + gge(7:18,7:18);
      H(rdof+[1:12])             = H(rdof+[1:12])             + hhe(7:18);
      K(rdof+[1:12])             = K(rdof+[1:12])             + kke(7:18);
      C(rdof+[1:12],rdof+[1:12]) = C(rdof+[1:12],rdof+[1:12]) + cce(7:18,7:18);
   end

end

v = x(Ndof+[1:Ndof]);

% Damping force.
D = C*v;

dxdt = zeros(Nx,1);
dxdt(1:Ndof) = v;
RHS = -G*v + H - K - D + F;
dxdt(Ndof+[1:Ndof]) = M\RHS;
dxdt(Nx) = WW;

fid = fopen('stat.txt','a');
fprintf(fid,'%+5.6e %+5.6e %+5.6e %+5.6e %+5.6e %+5.6e\n', ...
        t,x(Ndof-6+3),D(Ndof-6+1),D(Ndof-6+3),K(Ndof-6+1),K(Ndof-6+3));
fclose(fid);