function [Lmat,Rvec,y] = structureNL (x,t,s,P,F,shape,mdamp,grav, ...
                                      ret,slv,qg)
%
% Assemble the structural equations in state space form.  Note that
% q = (Lambda shape) eta.
%
%   States:           y vector:
%   eta       N       q            1:Ndj
%   deta/dt   N       dqdt     Ndj+1:2*Ndj
%                     d2qdt2 2*Ndj+1:3*Ndj
%
% Version:        Changes:
% --------        -------------
% 03.02.2018      Original code.
%
% Version:        Verification:
% --------        -------------
% 03.02.2018      NREL 5 MW modes match published values.
%
% Inputs:
% -------
% x               : States, eta, deta/dt, and qs.
% t               : time, not used.
% P               : Undeformed nodal positions.
% F               : Ndj nodal forces.
% shape           : Nret-by-Nmod mode shape matrix (= I if not used).
% freq            : Body mode frequencies, Hz, corresponding to the
%                   modes in shape.  Used for damping. C = 4*pi*zeta*M*f.
% grav            : Gravity vector in global coordinates.
% ret,slv         : List of retained and slave DOFs.
% qg              : A best guess for the slave DOFs.
%
% Outputs:
% --------
% Lmat, Rvec      : The rate of change of states dx/dt is Lmat\Rvec.
% y               : [q dq/dt d2q/dt2] vector.

%'structureNL'

N = size(shape,2);

[q,dqdt] = qFromx (x,s,P,shape,ret,slv,qg);
Nret = size(ret,1);
Nslv = size(slv,1);

% Build the nonlinear terms.
[M,MG,R,Q,Lambda,Gamma,Leq,dLdq] = buildMRQNL (s,q,dqdt,P,F);

gvec = zeros(Nret,1);
gvec(1:3) = grav;   % Foundation ref. node O_B/B0^g acceleration.
Qgrav = M*gvec;

Mr = (shape.')*M*shape;
RHS = (shape.')*(R + Q - MG + Qgrav);

if (real(s.zdamp) > 0)
   RHS = RHS - mdamp.*x(N+[1:N]);
end

sTs = (shape.')*shape;
Lmat = [sTs sparse(N,N);sparse(N,N) Mr];
Rvec = [sTs*x(N+[1:N]);RHS];

dxdt = Lmat\Rvec; % Needed for y output, but careful of conditioning.

% Accelerations.
dofs = [ret;slv];
Lamq = sparse(size(Lambda,1),size(Lambda,2));
Lamq(dofs,:) = Lambda;         % Undo the partitioning, so q = Lamq*qhat.
LSh = Lamq*shape;
Gamq = sparse(size(Gamma,1),size(Gamma,2));
Gamq(dofs,:) = Gamma;          % d2q/dt2 = Lamq*d2qh/dt2 + Gamq*dqh/dt.
Gsh = Gamq*shape;
d2qdt2 = LSh*dxdt(N+[1:N]) + Gsh*dxdt(1:N);

y = [q;dqdt;d2qdt2];


