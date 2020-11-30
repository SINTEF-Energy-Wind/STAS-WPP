%
% Newton-Raphson solution of the TC RWP grid with turbines, and
% generating linearized system matrices.
%
% S == 12 ==  8 ==  4 --  3 --  2 --  1           |  S  : substation node
%                                                 |  N  : turbine node
%               /-  7 --  6                       |  -- : collGrid1 cable
% S == 16 == 11                                   |  == : collGrid2 cable
%               \- 10 --  5
%
%         /- 14 --  9
% S == 15 
%         \- 18 -- 13 -- 17
%
%               /- 22 -- 21
% S == 20 == 19
%               \- 23 -- 26 -- 25
%  
%               /- 27 -- 30 -- 29
% S == 24 == 28
%               \- 32 -- 31
%
%
% S == Tp(66) == Ts(220) =cable= Tp(220) == PCC:Ts(440)
%

clear;

[length,time,mass,current,voltage,        ...
 velocity,force,power,stress,ndens,nvisc, ...
 stiffness,damping,resistance,inductance, ...
 capacitance,flux] = getLTMnorms ('LTMnorms.txt');

inpnm = '_P060_';
outnm = '_TCRWP_P060';

tnod = [1:32].';   % Let the node ID equal the turbine ID.
Nt   = size(tnod,1);
nsub = tnod(Nt) + 1;
npcc = nsub + 3;
Dia = 2*89.15;

Pe = 6e6*ones(Nt,1)                        / power;   % [MW], electrical power produced at each WT, here the same for all
we   = 50*2*pi                             * time;    % [rad/s], grid frequency
sq3  = sqrt(3);

% Nodes connected by cables
nod = [ 1  2; ...  % cable 1 connects node 1 to node 2
        2  3; ...
        3  4; ...
        4  8; ...
        5 10; ...
        6  7; ...
        7 11; ...
        8 12; ...
        9 14; ...
       10 11; ...
       11 16; ...
       12 nsub; ...
       13 18; ...
       14 15; ...
       15 nsub; ...
       16 nsub; ...
       17 13; ...
       18 15; ...
       19 20; ...
       20 nsub; ...
       21 22; ...
       22 19; ...
       23 19; ...
       24 nsub; ...
       25 26; ...
       26 23; ...
       27 28; ...
       28 24; ...
       29 30; ...
       30 27; ...
       31 32; ...
       32 28].';

Ncab = size(nod,2);
Nnod = max(max(nod));

load 'TurbinePosition.dat';  % x/D, y/D.
npos = [TurbinePosition;[17.5,15]]*Dia;
lencab = sqrt((npos(nod(2,:),1) - npos(nod(1,:),1)).^2 ...
       +      (npos(nod(2,:),2) - npos(nod(1,:),2)).^2);
typcab = [1 1 1 2 1 ...
          1 1 2 1 1 ...
          2 2 1 1 2 ...
          2 1 1 2 2 ...
          1 1 1 2 1 ...
          1 1 2 1 1 ...
          1 1].';

% Collection grid-66kV, type 1x3x95 Cu:
cableTypes.collGrid{1}.Ac   = 95;      % [mm^2], conductor area
cableTypes.collGrid{1}.L    = 0.440;   % [mH/km], inductance
cableTypes.collGrid{1}.C    = 0.170;   % [µF/km], capacitance
cableTypes.collGrid{1}.R_DC = 0.1933;  % [Ohm/km], AC resistance at 50 Hz
cableTypes.collGrid{1}.R_AC = 0.2466;  % [Ohm/km], DC resistance
% Collection grid-66kV, type 1x3x630 Cu:
cableTypes.collGrid{2}.Ac   = 630;     % [mm^2], conductor area
cableTypes.collGrid{2}.L    = 0.330;   % [mH/km], inductance
cableTypes.collGrid{2}.C    = 0.320;   % [µF/km], capacitance
cableTypes.collGrid{2}.R_DC = 0.0302;  % [Ohm/km], AC resistance at 50 Hz
cableTypes.collGrid{2}.R_AC = 0.0361;  % [Ohm/km], DC resistance
% Export sub. -220kV, type 1x3x1200 Cu:
cableTypes.exportSub.Ac     = 1200;    % [mm^2], conductor area
cableTypes.exportSub.L      = 0.361;   % [mH/km], inductance
cableTypes.exportSub.C      = 0.194;   % [µF/km], capacitance
cableTypes.exportSub.R_DC   = 0.0151;  % [Ohm/km], DC resistance at 50 Hz
cableTypes.exportSub.R_AC   = 0.0211;  % [Ohm/km], AC resistance
% Export land -220kV, type 3x1x1400 Cu:
cableTypes.exportLand.Ac    = 1400;    % [mm^2], conductor area
cableTypes.exportLand.L     = 0.510;   % [mH/km], inductance
cableTypes.exportLand.C     = 0.220;   % [µF/km], capacitance
cableTypes.exportLand.R_DC  = 0.0129;  % [Ohm/km], DC resistance at 50 Hz
cableTypes.exportLand.R_AC  = 0.0164;  % [Ohm/km], AC resistance

%- Offshore transformer (ot):
Pot  = 1.8e8*2                             / power;   % rating, 2 transformers
Votp = 66000                               / voltage; % voltage of primary winding
Vots = 220000                              / voltage; % voltage of secondary winding
Iotp = (Pot/3)/(Votp/sq3);
Lot  = 0.12*(Votp/sq3)/(Iotp*we);                     % inductance

%- Land transformer (lt):
Plt  = 1.8e8*2                             / power;   % rating, 2 transformers
Vltp = 220000                              / voltage; % voltage of primary winding
Vlts = 400000                              / voltage; % voltage of secondary winding
Iltp = (Plt/3)/(Vltp/sq3);
Llt  = 0.12*(Vltp/sq3)/(Iltp*we);

%- Export cable onshore (ons)
dons  = 30000                              / length;              % [m], export cable length onshore
Rlons = cableTypes.exportLand.R_AC*1e-3    / (resistance/length); % Export cable land-220kV
Llons = cableTypes.exportLand.L*1e-6       / (inductance/length);
Clons = cableTypes.exportLand.C*1e-9       / (capacitance/length);

%- Export cable offshore (off)
doff  = 30000                              / length;              % [m], export cable length offshore
Rloff = cableTypes.exportSub.R_AC*1e-3     / (resistance/length); % Export cable sub.-220kV
Lloff = cableTypes.exportSub.L*1e-6        / (inductance/length);
Cloff = cableTypes.exportSub.C*1e-9        / (capacitance/length);

%- Electrical export system
a12  = Votp/Vots;                                        % [-],   Offshore transformer Np/Ns turns ratio
I12  = [0 2000 4000 6000]                  / current;    % [A],   RMS phase currents for nonlinear inductance, for spline representation
L12  = [Lot Lot Lot Lot];                                % [H],   Lp + (a^2)Ls nonlinear phase inductances, for spline representation
R12  = 0.010                               / resistance; % [Ohm], Phase resistance - Check
L2   = 1                                   / inductance; % [H],   Effective inductance,  Bus 2 -> ground - Not used in TCRWP, but must be well-conditioned as 1/L2.
R2   = 0                                   / resistance; % [Ohm], Effective resistance,  Bus 2 -> ground - Not used in TCRWP, zero OK.
C2   = Cloff*doff;                                       % [F],   Effective capacitance, Bus 2 -> ground  
L23  = Lloff*doff + Llons*dons;                          % [H],   Phase inductance in cable
R23  = Rloff*doff + Rlons*dons;                          % [Ohm], Phase resistance in cable
L3   = 0.8132                              / inductance; % [H],   Effective inductance,  Bus 3 -> ground,  Shunt reactor
R3   = 0.510                               / resistance; % [Ohm], Effective resistance,  Bus 3 -> ground
C3   = Clons*dons;                                       % [F],   Effective capacitance, Bus 3 -> ground
a34  = Vltp/Vlts;                                        % [-],   Land transformer Np/Ns turns ratio
I34  = [0 2000 4000 6000]                  / current;    % [A],   RMS phase currents for nonlinear inductance, for spline representation
L34  = [Llt Llt Llt Llt];                                % [H],   Lp + (a^2)Ls nonlinear phase inductances, for spline representation
R34  = 0.010                               / resistance; % [Ohm], Phase resistance - Check

Rlen = zeros(Ncab,1);
Llen = zeros(Ncab,1);
Clen = zeros(Ncab,1);
for idx = 1:Ncab
  cabTypeNum = typcab(idx);
  Rlen(idx)  = cableTypes.collGrid{1,cabTypeNum}.R_AC*(1e-3) / (resistance/length);
  Llen(idx)  = cableTypes.collGrid{1,cabTypeNum}.L   *(1e-6) / (inductance/length);
  Clen(idx)  = cableTypes.collGrid{1,cabTypeNum}.C   *(1e-9) / (capacitance/length);
end

rand('state',1);
lencab = lencab.*(0.99 + 0.02*rand(size(lencab)));  % Prevent perfect symmetry.

Rcab = Rlen.*lencab; %resistance
Lcab = Llen.*lencab; %inductance

Cnod = zeros(Nnod,1);
Cnod(1:Nnod-1) = 0.5*Clen;
Cnod(2:Nnod) = Cnod(2:Nnod) + 0.5*Clen;

Nu = 3 + 2*Nnod;  % u = [we, vpccd,q, itd,q for all nodes]
Nxf = 2*Ncab+4*Nnod+14;
gpar = [Rcab;Lcab;Cnod];
xpar = [a12 I12 L12 R12 L2 R2 C2 L23 R23 L3 R3 C3 a34 I34 L34 R34].';

ldofs = [2*Ncab+2*Nnod+[1:2*Nnod] 2*Ncab+4*Nnod+[3:4]].';
[ip,ret,cr] = partitionMatrix ([1:Nxf].',ldofs,[]);
Ndel = size(ldofs,1);
Nret = size(ret,1);

% Solve the load flow problem.  Solve the combined equations dx/dt = 0, and
% P = id*vd+iq*vq, Q = id*vq - iq*vd or
% P = |i||v| cos(thv-thi), Q = |i||v| sin(thv-thi)
Phat = [Pe; 0];
Qhat = [zeros(Nt,1); 0];

i2na = [1:2:2*Nnod-1].';
i2nb = [2:2:2*Nnod].';
i2ca = [1:2:2*Ncab-1].';
i2cb = [2:2:2*Ncab].';

Inom = Phat./Votp;  %nominal current at nodes

u0 = zeros(Nu,1);
u0(1:3) = [we;Vlts;0]; %secondary terminal of land tranformer is set to voltage of regional grid v4d,q
u0(3+i2na) = Inom;  % Nominal zero-offset power inflow.

xf0 = zeros(Nxf,1);

intNod = setdiff(nod(1,:),nod(2,:)); %find initial nodes that are end points of the grid
xf0(2*intNod-1) = Inom(intNod); %current in cables connecting initial nodes to the rest of the grid

missingNodes  = setdiff([1:nsub],[intNod, nsub]); %missing, i.e. not connected nodes (excluding initial nodes and substation node)
conNod = [];
for jj = intNod
  idxNod = jj;   
    while sum([nod(2,:) == nod(2,idxNod)]) == 1  &&  ~isempty(missingNodes)
      nextNod = nod(2,idxNod); %node connected to node jj
      xf0(2*nextNod-1) = xf0(2*idxNod-1) + Inom(nextNod);
      idxNod = nod(2,idxNod);
      missingNodes(find(missingNodes == idxNod)) = [];
    endwhile
    nextNod = nod(2,idxNod); %node connected to node jj
    conNod = [conNod; nextNod];
end

startNod = unique(setdiff(conNod,nsub));  %use not yet connected nodes as start, except nsub
for jj = startNod'
  if isempty(intersect(missingNodes,nod(1,find([nod(2,:) == jj])))) %ensure that nodes connected from grid endpoints to idxNod have been included
    idxNod = nod(1,find([nod(2,:) == jj]));
    nextNod = jj; %node connected to node jj
    xf0(2*nextNod-1) = sum(xf0(2*idxNod-1)) + Inom(nextNod);
    missingNodes(find(missingNodes == jj)) = [];  %was included above
    
    if ~isempty(missingNodes)
      %Continue as above with cables without branches
      idxNod = jj; %start from connection node
      while sum([nod(2,:) == nod(2,idxNod)]) == 1 
        nextNod = nod(2,idxNod); %node connected to node jj
        xf0(2*nextNod-1) = xf0(2*idxNod-1) + Inom(nextNod);
        idxNod = nod(2,idxNod);
        missingNodes(find(missingNodes == idxNod)) = [];
      endwhile
      nextNod = nod(2,idxNod); %node connected to node jj
      conNod = [conNod; nextNod];
      end

  else
    warning(['Not all cables from endpoints of the grid that branch from ',num2str(jj),' have been handled.'])
  end
  
end

if ~isempty(missingNodes)
  warning(['Not all nodes have been connected to the grid. Missing nodes: ',num2str(missingNodes)])
else
  clear missingNodes
end

xf0(2*Ncab+i2na) = Votp;  % Voltage at grid nodes

%--------------------------------- %
% ------- AC export system ------- %
%--------------------------------- %
%   i12d,q   1,2   - offshore transformer
%   ir2d,q   3,4   - offshore compensation
%   v2d,q    5,6   - offshore transformer, secondary voltage, at side of export cable
%   i23d,q   7,8   - export cable
%   ir3d,q   9,10  - shunt reactor
%   v3d,q   11,12  - land transformer, primary voltage, at side of export cable
%   i34d,q  13,14  - land transformer
i12d0 = sum(xf0(find(nod(2,:)==nsub)*2-1));   %initial i12d as sum of currents in cables to node nsub
xf0(2*Ncab+4*Nnod+[1:14]) = [i12d0;0; ...      % i12d,q is Inom in cable nsub (from node nsub)
                             0;0; ...
                             Vots;0; ...
                             a12*i12d0;0; ...  % current in export cable; i12d,q corrected by offshore transformer turns ratio
                             0;0; ...
                             Vltp;0; ...
                             a12*i12d0;0];     % shunt reactor compensates 100% of both cables' reactive power 
                             
                             
x0 = xf0(ret);

%{
del = sqrt(eps);
[dx0dt,A,B] = buildWPPgrid (1,x0,u0,gpar,xpar,nod,nsub);
Ac = zeros(size(A,1),size(A,2));
Bc = zeros(size(A,1),size(u0,1));
for ix = 1:size(A,1)
   xc = x0;
   xc(ix) = xc(ix) + i*del;
   [dxcdt,jnkA,jnkB] = buildWPPgrid (0,xc,u0,gpar,xpar,nod,nsub);
   Ac(:,ix) = imag(dxcdt)/del;
end
sparse(A - Ac)
for iu = 1:size(B,2)
   uc = u0;
   uc(iu) = uc(iu) + i*del;
   [dxcdt,jnkA,jnkB] = buildWPPgrid (0,x0,uc,gpar,xpar,nod,nsub);
   Bc(:,iu) = imag(dxcdt)/del;
end
sparse(B-Bc)
return
%}

[dx0dt,A,B] = buildWPPgrid (0,x0,u0,gpar,xpar,nod,nsub);  %nonlinear equations
iuied = 3 + i2na; %itd
iuieq = 3 + i2nb; %itq
ixvd  = 2*Ncab + i2na;  %vds
ixvq  = 2*Ncab + i2nb;  %vq
P = u0(iuied).*x0(ixvd) + u0(iuieq).*x0(ixvq);
Q = u0(iuied).*x0(ixvq) - u0(iuieq).*x0(ixvd);

Res = [dx0dt;P-Phat;Q-Qhat];
Rval = (Res.')*Res;

z = [x0;u0(3+[1:2*Nnod])];  % [iLd,q; vd,q; i12d,q; v2d,q; i23d,q; ir3d,q; v3d,q; i34d,q; ... % x
                            %  itd,q]                                                         % u
Nz = size(z,1);

% ===================================================================
% Newton-Raphson method solution.
cnv = eps^0.4;
Ns = 100;
%bta = [0.01;0.02;0.05;0.1;0.2;0.5;ones(Ns,1)];
bta = ones(Ns,1);
litmax = 20;
conv = 0;
iter = 0;
Niter = 0;
lam = 1;
while (((real(Rval) > cnv) && (iter < Ns)) || (iter == 0))
   iter = iter + 1;

   % Compute the tangent function at the latest point.
   x = z(1:Nret);
   u = [u0(1:3);z(Nret+[1:2*Nnod])];  % update inputs itd,q with states in z
   [dxdt,A,B] = buildWPPgrid (1,x,u,gpar,xpar,nod,nsub);  %linear equations

   dRdz = zeros(Nz,Nz);
   
   ir = [1:Nret]; %x
   ic = ir;
   dRdz(ir,ic) = A;
   ic = Nret+[1:2*Nnod];
   dRdz(ir,ic) = B(:,3+[1:2*Nnod]);
   
   ir = Nret+[1:Nnod]; %P
   ic = ixvd;
   dRdz(ir,ic) = diag(u(iuied));
   ic = ixvq;
   dRdz(ir,ic) = diag(u(iuieq));
   ic = Nret + i2na;
   dRdz(ir,ic) = diag(x(ixvd));
   ic = Nret + i2nb;
   dRdz(ir,ic) = diag(x(ixvq));
   
   ir = Nret+Nnod+[1:Nnod]; %Q
   ic = ixvd;
   dRdz(ir,ic) = -diag(u(iuieq));
   ic = ixvq;
   dRdz(ir,ic) = diag(u(iuied));
   ic = Nret + i2na;
   dRdz(ir,ic) = diag(x(ixvq));
   ic = Nret + i2nb;
   dRdz(ir,ic) = -diag(x(ixvd));

   printf('   Matrix condition %+6.3e\n',cond(dRdz));
   fflush(stdout);

%[slap,shp,ifrq] = eigVal (dRdz);

%return

   dz = -dRdz\Res;  % -(df/dz)^-1 * f(z)
   lflg = 0;
   liter = 0;
   while ((lflg == 0) && (liter < litmax))
      liter = liter + 1;

      z1 = z + bta(iter)*lam*dz; % z1 = z + Dz, with modified Newton-Raphson iteration Dz = - alpha * (df/dz)^-1 * f(z)

      x = z1(1:Nret);
      u = [u0(1:3);z1(Nret+[1:2*Nnod])];    % update inputs itd,q with states in z
      [dx1dt,jnkA,jnkB] = buildWPPgrid (0,x,u,gpar,xpar,nod,nsub);  %nonlinear equations
      P1 = u(iuied).*x(ixvd) + u(iuieq).*x(ixvq);
      Q1 = u(iuied).*x(ixvq) - u(iuieq).*x(ixvd);

      Res1 = [dx1dt;P1-Phat;Q1-Qhat];  % f(z)
      R1 = (Res1.')*Res1;              % Residual R := f^T f

      printf('   %5d %5d  %+5.3e  %+5.3e\n',iter,liter,Rval,R1);
      fflush(stdout);

      if (real(R1) < Rval)
         % OK!  Prepare for the next iteration.
         lflg = 1;
         Res = Res1;
         Rval = R1;
         z = z1;
      else
         % Backtrack.
         lam = 0.5*lam;
         if (liter == litmax)
            [iter R1 Rval]
            printf('Warning, proceeding without lambda convergence.\n');
            lflg = 1;
            Res = Res1;
            Rval = R1;
            z = z1;
return
         end
      end

   end  % Newton inner.

   if (iter == Ns) && (real(Rval) > cnv)
      Rval
      printf('Warning, max iterations, proceeding without Rval convergence.\n');
return
   end

end

xf = zeros(Nxf,1);
xf(ret) = x;

Ag = A;
Bg = B;
xgf = xf;
Bgq = -Ag\Bg;
gret = ret;
eval(["save('-binary','Ag" outnm ".bin','Ag');"]);
eval(["save('-binary','Bg" outnm ".bin','Bg');"]);
eval(["save('-binary','xgf" outnm ".bin','xgf');"]);
eval(["save('-binary','Bgq" outnm ".bin','Bgq');"]);
eval(["save('-binary','gret" outnm ".bin','gret');"]);


