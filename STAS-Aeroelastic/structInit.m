function [xs,etas,q,dq,d2q,P,shape,sh2,freq,mdamp,ret,slv, ...
          Ndj,Nnod,Neta] = structInit (s,yaw,psi,betas,Omega,modeflag)
%
% Initialize parameters related to the structural model.
%
% Version:        Changes:
% --------        -------------
% 25.04.2018      Original code.
%
% Version:        Verification:
% --------        -------------
% 25.04.2018      
%
% Inputs:
% -------
% s               : Data structure.
% yaw,psi,betas   : Yaw, azimuth, blade pitch angles.
% Omega           : Rotor speed.
% modeflag        : 1 for structural modes, 0 for full DOFs.
%
% Outputs:
% --------
% xs              : Constrained [qhat;dqhat/dt]
% etas            : Reduced [eta;deta/dt]
% q,dq,d2q        : Structural DOF displacement, velocity, accel.
% P               : Nodal offsets.
% shape,freq      : Mode shapes and frequencies.
% mdamp           : Vector of modal damping ratios.
% ret,slv         : Retained and slave DOFs.
%

%'structInit'

Nnod = s.foundation.Nnod + s.tower.Nnod + s.nacelle.Nnod     ...
     + s.driveshaft.Nnod + s.blade(1).Nnod + s.blade(2).Nnod ...
     + s.blade(3).Nnod;

[idofs,idofm,inods,inodm,Ndof] = getDOFRefs (s);
[Tn_y,Th_d,Tb_h] = basicTransforms (s.nacelle.delta,s.driveshaft.phi);
Pin = assemblePin (s);
[q,P,Ts0_B,TB0_g] =                                                      ...
      undeformedPosition (Pin,yaw,s.nacelle.delta,psi,s.driveshaft.phi, ...
                          betas,0,idofs,idofm,inods,inodm);

% Initial velocity considering rotor rotation; this fills in the 
% state vector for structural DOFs, xs, including both positions and
% velocities.  xs includes the nominal constraint equations, but not
% locked joints or modal DOF reduction.
[xs,dq,ret,slv] = initialDOFVelocity (s,q,P,Omega);

% Ndj: the total number of structural DOFs, including the six joint
% DOFs.
Ndj  = size(q,1);
Nret = size(ret,1);
Nslv = size(slv,1);

% Initial accelerations can be obtained via structureNL outputs.
% We don't yet know the applied forces, let these be zero.  The 
% output of structureNL is the initial y vector for the structural
% DOFs.
%
% [Caution, this won't work properly with s.zdamp; not an issue
%  here if we don't use d2q.]

Force = zeros(Ndj,1);
%[Lms,Rvs,ys] = structureNL (xs,0,s,P,Force,speye(Nret),sparse(Nret,1), ...
%                            ret,slv,q(slv));
d2q = zeros (Ndj,1); % ys(2*Ndj+[1:Ndj]);  % Zero acceleration for
                                           % computing body modes.

% And get the mode shapes.
[M,dM,MG,dMG,dMGd,R,dR,dRd,Q,dQ,dQd,slv,ret, ...
 Lambda,Gamma,Leq,dLdq,dL,dGam,dGamd] =      ...
                    buildMRQLin (s,q,dq,d2q,P,Force);
K = dM + dMG - dR - dQ;

%dM
%dMG
%dR
%dQ

if (modeflag == 1)
   [shape,freq,mdamp] = bodyModes (s,M,K,ret,slv);
   Neta = size(shape,2);
elseif (modeflag == 0)
   shape = speye (Nret);
   freq  = zeros (Nret,1);
   mdamp = zeros (Nret,1);
   Neta  = Nret;
end
sh2 = [shape sparse(Nret,Neta);sparse(Nret,Neta) shape];

etas = (sh2.')*xs;

