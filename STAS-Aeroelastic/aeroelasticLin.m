function [LHS,A,Bu,By,C,Du,Dy,shape,freq,mdamp,ret,slv,Psi] =     ...
                  aeroelasticLin (x,Vg,s,a,q0,dq0dt,d2q0dt2,P,F0, ...
                                  psiFlag,modflag,shpflag,        ...
                                  shapein,mdampin,grav)
%
% Build a linear model of the turbine at a given operating point, which
% does not need to be an equilibrium point.  In the output matrices,
% joints are free, and need to be controlled, or locked by eliminating
% their DOFs.
%
%   States:              y vector:             u vector:
% ----------------------- Structure --------------------------
%   eta        N         q         Ndj         F          Ndj
%   deta/dt    N         dq/dt     Ndj
%                        d2q/dt2   Ndj
%                        xng     3*Nnod
%                        vng     6*Nnod
%                        F         Ndj
% ---------------------- Aerodynamic -------------------------
%   ad         1         (q)       1:24        Vg          3
%   a1,a2     2:3        (dq/dt)  25:48
%   Vih z,t   4:5        xng1,2   49:54
%   Vi z,t    6:7        vng1,2   55:60
%                        wg       61:63
%                        wa       64:66
%                        aq        67
%                        Cl,Cd,Cm 68:70
%                        Fl,Fd,M  71:73
%                        Fa       74:79
%                        Fp       80:85
%                        Fr       86:91
%                        Fzts     92:94
%                        Vg       95:97
%                        Ua       98:100
%                        Umag      101
%                        Ur      102:104
%                        Uzts    105:107
%                        Vzts    108:110
%                        Viq     111:112
%                        Viy     113:114
%                        Vixyz   115:117
%                        Wmag      118
%                        xeg     119:121
%                        xhg     122:124
%                        xnr1,2  125:130  
%                        xer     131:133
%                        r         134
%                        Lp        135
%                        z         136
%                        f         137  
% (Repeat the above x,y,u consecutively for each blade element.)
%                        Azi       (1)  (effective rotor azimuth)
%                        Waero     (1)  (aero rotor speed)
%                        Dp        (3)
%
%
% Version:        Changes:
% --------        -------------
% 06.02.2018      Original code.
%
% Version:        Verification:
% --------        -------------
% 06.02.2018      
%
% Inputs:
% -------
% ns              : Number of structural, aerodynamic, electrical, and
%                   control states.
%
% Outputs:
% --------
% 

%'---------------------'
%'Lin'

Ndj = size(q0,1);

[idofs,idofm,inods,inodm,Ndof] = getDOFRefs (s);
Nnod = Ndof/6;
[Tn_y,Th_d,Tb_h] = basicTransforms (s.nacelle.delta,s.driveshaft.phi);

% =========== Structural state equations =============

% Add generator torque to the mean force vector.
[lls,aas,bbus,bbys,ccs,ddus,ddys,                     ...
 Lambda,Gamma,shape,freq,mdamp,ret,slv] =             ...
               structureLin (modflag,shpflag,         ...
                             s,q0,dq0dt,d2q0dt2,P,F0, ...
                             shapein,mdampin,grav);
Nxs  = size(aas,1);
Nus  = size(bbus,2);
Nysq = size(ccs,1);
Nys  = Nysq + 9*Nnod + Ndj;
Nxs2 = Nxs/2;

% =========== Aerodynamic state equations ============ 
if (a.icp(1) == 0)
   Ncp = a.Neb;
else
   Ncp = size(a.icp,1);
end
Nxa = 7*Ncp*a.Nb;
xaero = x(Nxs+[1:Nxa]);
bsh = bladeModeShape (s,ret,shape);
[dxadt,aaa,bbya,cca,ddya,Psi] = ...
          aeroLin (psiFlag,xaero,Vg,s,a,q0,dq0dt,P,bsh);
Nxa = size(aaa,1);
Nya = size(cca,1);
Nua = 3*a.Nb*a.Neb;  % We will define Vg as a global input.
Nyae = 137;          % Number of aero y's associated with each element.
Nea = a.Nb*a.Neb;
Neb = a.Neb;
lla = speye (Nxa);

% ======= Equations for linking aero-structure ======= 
[xng,vng,ddypv] = globalPosVel (s,idofs,q0(1:Ndof),P(1:Ndof),dq0dt(1:Ndof));
[Waero,azi,ddyRSA] = rotorSpeedAero (q0,dq0dt,P,Tn_y, ...
                                     idofs(3),idofs(4),idofm(6)-6);

[llas,aaas,bbuas,bbyas,ccas,dduas,ddyas] =                       ...
               linkAeroelastic (s,a,                             ...
                                lls,aas,bbus,bbys,ccs,ddus,ddys, ...
                                lla,aaa,bbya,cca,ddya,           ...
                                ddypv,ddyRSA);

% (Temporary.)
LHS = llas;
  A = aaas;
 Bu = bbuas;

Nx = size(aaas,1);
Ny = size(ccas,1);
Nu = size(Bu,2);
By = sparse(Nx,Ny);
C  = sparse(Ny,Nx);
Du = sparse(Ny,Nu);
Dy = sparse(Ny,Ny);

By(:,1:Ny) = bbyas;
C(1:Ny,:) = ccas;
Du(1:Ny,:) = dduas;
Dy(1:Ny,1:Ny) = ddyas;
