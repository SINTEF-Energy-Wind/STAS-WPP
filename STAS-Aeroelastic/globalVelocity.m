function [vg0,Dy] = globalVelocity (qn,qB,Pn,PB,dqndt,dqBdt)
%
% Computes the nodal velocity [v;w] and the Dy matrix entries for linear
% perturbations.
%
%   States:           y vector:         u vector:
%                     qB
%                     qn
%                     dqBdt
%                     dqndt
%                     [vg,wg]
%
% Version:        Changes:
% --------        -------------
% 20.12.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 20.12.2017      Verified using complex step via vg0 output.
%
% Inputs:
% -------
% qn              : DOFs associated with the node.
% qB              : DOFs associated with the body reference node.
% Pn              : Undeformed pos,rot of the node wrt body reference.
% PB              : Undeformed pos,rot of the body reference.
% dqndt,dqBdt     : Rates of change of qn, qB.
%
% Outputs:
% --------
% vg0             : Nodal velocities at the input conditions.
% Dy              : dv/dq and dv/dqdot, 6-by-24 for each node.

%id = tic();

Nnod = size (qn,1)/6;

vg0 = zeros(6*Nnod,1);
Dy  = zeros(6,24*Nnod);

for inod = 1:Nnod

%inod
%toc(id)

   jc24 = 24*(inod-1);
   jc12 = 12*(inod-1);
   jc6 = 6*(inod-1);
   Qg  = Qgnod (qn,qB,Pn,PB);
%'Qg'
%toc(id)
   dQg = dQgdq (qn,qB,Pn,PB);
%'dQg'
%toc(id)
   vg0 = Qg*[dqBdt(jc6+[1:6]);dqndt(jc6+[1:6])];

   Dy(:,jc24+[13:24]) = Qg;
   for idof = 1:6
      kc12 = 12*(idof-1);
      Dy(:,jc24+idof) = dQg(:,kc12+[1:12])*[dqBdt(jc6+[1:6]);dqndt(jc6+[1:6])];
   end
   for idof = 1:6
      kc12 = 12*(idof+6-1);
      Dy(:,jc24+6+idof) = dQg(:,kc12+[1:12])*[dqBdt(jc6+[1:6]);dqndt(jc6+[1:6])];
   end
%toc(id)
end

%{
'---globalVelocity---'
[qn qB Pn PB dqndt dqBdt]
Qg
dqBdt(jc6+[1:6])
dqndt(jc6+[1:6])
vg0
'--end globalVelocity--'
%}

%'end globalVelocity'
%toc(id)