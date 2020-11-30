function [xs,etas,q,qd,qdd,P,shape,freq,mdamp,ret,slv, ...
          Ndj,Nnod,Neta] = structInit (s,yaw,psi,betas,Omega)
%
% Initialize parameters related to the structural model.
%
% The mode shape calculation should ideally include centrifugal
% stiffening, but this requires solving for the displaced position
% of the blade under the rotational acceleration field.
% - start the rotor spinning
% - isolate one of the blades, fix its root
% - compute the nonlinear steady-state solution
% - compute the eigenmodes.
% This is left to a future implementation.  Note that the primary
% purpose of the modal reduction is to prevent very high-frequency
% structural modes from making system matrices ill-conditioned.
% A modal reduction based on static-blade modes, although not
% ideal, is sufficient for this purpose.
%
% Version:        Changes:
% --------        -------------
% 25.04.2018      Original code.
% 25.01.2020      Accommodate updated buildMRQ.
%
% Version:        Verification:
% --------        -------------
% 25.04.2018      
% 25.01.2020      
%
% Inputs:
% -------
% s               : Data structure.
% yaw,psi,betas   : Yaw, azimuth, blade pitch angles that 
% Omega           : Rotor speed.
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

Nnod = s.foundation.Nnod + s.tower.Nnod + s.nacelle.Nnod     ...
     + s.driveshaft.Nnod + s.blade(1).Nnod + s.blade(2).Nnod ...
     + s.blade(3).Nnod;

[idofs,idofm,inods,inodm,Ndof] = getDOFRefs (s);
[Tn_y,Th_d,Tb_h] = basicTransforms (s.nacelle.delta,s.driveshaft.phi);
Pin = assemblePin (s);
[q,P,Ts0_B,TB0_g] =                                                      ...
      undeformedPosition (Pin,yaw,s.nacelle.delta,psi,s.driveshaft.phi, ...
                          betas,0,idofs,idofm,inods,inodm);
Ndj  = size(q,1);

% Initial velocity considering rotor rotation.
[xs,qd,ret,slv] = initialDOFVelocity (s,q,P,Omega);
Nret = size(ret,1);

% We don't yet know the applied forces, let these be zero.  The
% DOF accelerations should be set to zero for the initial body mode
% calculation.  (After a static solution to the state equations is
% available the body modes may be recomputed so as to include
% centrifugal stiffening and other effects.)
F     = zeros(Ndj,1);
grav  = zeros(3,1);
qdd   = zeros (Ndj,1); % Zero acceleration for computing body modes.

% Get the unconstrained matrices.
[M,G,H,C,K,Q,dM,dMg,dG,dGd,dH,dHd,dC,dCd,dK,dQ,dQd] = ...
                     buildMRQ (1,s,q,qd,qdd,P,F,grav);
KK = dG - dH + dK;  % Should dG and dH be included, since the static
                    % solution has not been obtained?

% Constrain.
qh   = xs(1:Nret);
qhd  = xs(Nret+[1:Nret]);
qhdd = qdd(ret);
[q,qd,qdd,Tqhq,dTdq,dTdqd] = buildTqhq (0,s,qh,qhd,qhdd,P,ret,slv);
T = Tqhq(Ndj+[1:Ndj],Nret+[1:Nret]);
Mr   = (T.')*M*T;
Kr   = (T.')*KK*T;

[shape,freq,mdamp] = bodyModes (s,Mr,Kr,ret,slv);
Neta = size(shape,2);

sh2 = [shape sparse(Nret,Neta);sparse(Nret,Neta) shape];
etas = (sh2.')*xs;  % Best estimate.
