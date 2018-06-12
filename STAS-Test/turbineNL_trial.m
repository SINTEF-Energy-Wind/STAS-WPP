function [Lmat,Rvec,y] = turbineNL_trial (t,x1,y1,            ...
                                    psiFlag,s,a,P,            ...
                                    ret,slv,shape,mdamp,grav, ...
                                    Tas,Try,Vg,ch,Lel,foilwt, ...
                                    aoaz,aoast,xas,yas,Psi)
%
% Build the nonlinear turbine state equations (giving an overall dx/dt
% and updated output y) at a specified point in time, where a state
% and output vector are available.
%
%   States:           y vector:
%   eta       N       q            Ndj
%   deta/dt   N       dqdt         Ndj
%                     d2qdt2       Ndj
%
% (For each aero element)
%   ad         1      Fa            6
%   a1,a2     2:3      
%   Vih z,t   4:5      
%   Vi z,t    6:7
%
%
%                     F            Ndj
%                     Tgen          1
%
% Version:        Changes:
% --------        -------------
% 17.04.2018      Original code.
%
% Version:        Verification:
% --------        -------------
% 17.04.2018      
%
% Inputs:
% -------
% s,a             : Structural and aero data structures.
%
% Outputs:
% --------
% 

%'---------------------'
%'NL'

load 'LTMnorms.txt';
length  = LTMnorms(1);
time    = LTMnorms(2);
mass    = LTMnorms(3);
current = LTMnorms(4);
power   = mass*(length^2)/(time^3);
torque  = power*time;

Nxs = 2*size(shape,2);
Nxa = size(Psi,2);
Nx  = size(x1,1);
Ny  = size(y1,1);
Ndj = size(P,1);

[idofs,idofm,inods,inodm,Ndof] = getDOFRefs (s);

xs = x1(1:Nxs);
xa = x1(Nxs+[1:Nxa]);

q    = y1(1:Ndj);
dqdt = y1(Ndj+[1:Ndj]);

% Aerodynamic state equations.
[dxadt,ya] = aeroNL (psiFlag,s,a,xa,t,q,dqdt,P, ...
                     Tas,Try,Vg,ch,Lel,foilwt,  ...
                     aoaz,aoast,xas,yas,Psi);

% Control state equations.
Wmeas = dqdt(Ndj-4);

% Tjaereborg.
%Tgen = (((Wmeas - 2.2965)/(2.3415 - 2.2965))*2.2e6/2.3)/torque;

% NREL 5 MW.
Wtarget = 0.95;
Tgen = 5e7*(Wmeas - Wtarget)/torque;

% Update force vector.  The torque is applied at the rear bearing,
% and is defined about the local driveshaft master node orientation.
Fa = ya;
F = airfoilToNodal (s,a,Fa,Tas,q,P,idofs);

% Torque aligned with element orientation.
%qn1 = q(idofm(4)-6+[1:6]);
%qn2 = q(idofm(4)+[1:6]);
%Pn1 = P(idofm(4)-6+[1:6]);
%Pn2 = P(idofm(4)+[1:6]);
%[xenac,TsBnac] = elementCSFromNodes (qn1,qn2,Pn1,Pn2);
%F(idofs(4)+6) = F(idofs(4)+6) - Tgen;
%F(idofm(4)+[4:6]) = F(idofm(4)+[4:6]) + TsBnac*[-Tgen;0;0];

% Torque aligned with master node orientation.  
%
% Option 1, put the torque on the second node on the driveshaft,
% in order to avoid issues with how to properly apply an external
% force at the reference node.
Tnn0 = TFromTheta (q(idofm(4)+[4:6]));
Tn0B = TFromTheta (P(idofm(4)+[4:6]));
F(idofs(4)+12) = F(idofs(4)+12) - Tgen;
F(idofm(4)+[4:6]) = F(idofm(4)+[4:6]) + Tn0B*Tnn0*[-Tgen;0;0];

% Option 2, apply the torque across the joint DOF.
%F(Ndj-4) = F(Ndj-4) - Tgen;  % Check +/-.

% Structural state equations.
[Lms,Rvs,ys] = structureNL (xs,t,s,P,F,shape,mdamp,grav, ...
                            ret,slv,q(slv));

% Electrical state equations.
%dxedt = electricNL ();

% (Add additional systems here.)


Lmat = [Lms sparse(Nxs,Nxa);sparse(Nxa,Nxs) speye(Nxa)];
Rvec = [Rvs;dxadt];
y    = [ys;ya;F;Tgen];  % Now up to date with the input x.

