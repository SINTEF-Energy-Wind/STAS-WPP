function [M,C,K,slv,ret,Lambda,Gamma,                      ...
          Leq,dLdq,dLam,dGam,dGamd,shape,freq,mdamp] =     ...
                        buildMCK (shpflag,s,q,dqdt,d2qdt2, ...
                                  P,F,shapein,mdampin)
%
% Assemble the linearized mass, gyroscopic/damping, and stiffness
% matrices for the turbine structure.
%
% Version:        Changes:
% --------        -------------
% 23.11.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 23.11.2017      NREL 5 MW modes match published values.
%
% Inputs:
% -------
% shpflag         : 1 to compute shape, damping etc from the present
%                   matrices, 0 to use existing basis functions, via
%                   shapein and mdampin.
% s               : Wind turbine data structure.
% q...d2qdt2      : Initial displacements, velocities, accelerations.
% P               : Nodal positions.
% F               : Initial forces.
%
% Outputs:
% --------
% M,C,K           : Mass, damping, stiffness matrices.
% slv,ret         : Lists of slave and retained DOFs.
% Lambda          : DOF constraint matrix.
% Gamma           : Constraint force matrix.
% shape           : Body mode shapes.
% freq            : Body mode natural frequencies.

[M1,dM1,MG1,dMG1,dMGd1,R1,dR1,dRd1,Q1,dQ1,dQd1,slv,ret, ...
 Lambda,Gamma,Leq,dLdq,dLam,dGam,dGamd] =               ...
                      buildMRQLin (s,q,dqdt,d2qdt2,P,F);

C1 = dMGd1 - dRd1 - dQd1;
K1 = dM1 + dMG1 - dR1 - dQ1;

if (shpflag == 1)

   % Compute mode shapes for each body, based on the relevant portions
   % of the mass and stiffness matrices.
   [shape,freq,mdamp] = bodyModes (s,M1,K1,ret);

elseif (shpflag == 0)

   % Use existing basis functions and modal damping.
   shape = shapein;
   mdamp = mdampin;
   freq  = zeros(size(shape,2));

end

Nmdof = size(shape,2);

% Convert the matrices to modal DOFs.
M = (shape.')*M1*shape;
C = (shape.')*C1*shape;
K = (shape.')*K1*shape;

edofs = [1:Nmdof].';
damp = sparse (edofs,edofs,mdamp,size(C,1),size(C,2));
C = C + damp;
