function [q,dqdt] = qFromx (x,s,P,shape,ret,slv,qg)
%
% Compute the full nodal DOF vector and its rate of change based on
% the state input.
%
%   States:           y vector:
%   eta       N       q         1:Ndj
%   deta/dt   N       dqdt  Ndj+1:2*Ndj
%
% Version:        Changes:
% --------        -------------
% 16.02.2018      Original code.
%
% Version:        Verification:
% --------        -------------
% 16.02.2018      Checked some sample cases.  Checked that constraint
%                 equations are satisfied to high precision.
%
% Inputs:
% -------
% x               : States, eta, deta/dt.
% P               : Undeformed nodal positions.
% shape           : There are two options for the shape matrix.  One,
%                   it can be the body modes, an Nret-by-Nmod matrix.
%                   Two, it can be a global state-space mode shape
%                   matrix, 2*Nret-by-Nmod.  Use speye(Nret) if there
%                   is no modal reduction.
% ret,slv         : List of retained and slave DOFs.
% qg              : Initial estimate of qs.
%
% Outputs:
% --------
% q,dqdt          : Nodal DOFs.

%'qFromx'

N    = size(shape,2);
Nsh  = size(shape,1);
Nret = size(ret,1);
Nslv = size(slv,1);
Ndj  = Nret + Nslv;

[idofs,idofm,inods,inodm,Ndof] = getDOFRefs (s);
[Tn_y,Th_d,Tb_h] = basicTransforms (s.nacelle.delta,s.driveshaft.phi);

q    = zeros(Ndj,1);
dqdt = zeros(Ndj,1);

if (Nsh == Nret)
   % shape applies independently to position and velocity components
   % of x, to give retained q, dqdt DOFs.
   q(ret) = shape*x(1:N);
   dqdt(ret) = shape*x(N+[1:N]);
elseif (Nsh == 2*Nret)
   % shape multiplies x to give retained [q;dqdt] DOFs.
   qdq = shape*x;
   q(ret) = qdq(1:Nret);
   dqdt(ret) = qdq(Nret+1:2*Nret);
end
qs = qg;

% Let the constraint equations be written C(q) = 0, with first variation
% L(q) dq = 0.  The task is to find qs such that C(q) = 0 is satisfied.
% This can be accomplished by Newton's method, -L(q0) dq = C(q0).
% [Should add some robustness via the Press et al. algorithm.]
bta   = 1;
cnv   = eps^(0.75);  % Has to be below complex step precision, otherwise
Ns    = 100;         % no iteration occurs in the complex step case...
conv  = 0;           % [Perhaps fixed with iter <= 1 criterion.]
iter  = 0;
Niter = 0;
Rval  = 1;
while ((iter <= 1) || ((real(Rval) > cnv) && (iter <= Ns)))
   iter   = iter + 1;
   q(slv) = qs;
   [Lambda,L,C,jnkret,jnkslv] = constraints (q,P,Tb_h,idofs,idofm);
   Ls     = L(:,slv);
   dqs    = -Ls\C;
   qs     = qs + bta*dqs;
%   Rval   = max(abs(real(C)));  % Does not work with complex step.
   Rval   = max(abs(C));  % Works with complex step.
end

% Assemble dq/dt slave DOFs.
mLsL = Lambda(Nret+[1:Nslv],:);
dqdt(slv) = mLsL*dqdt(ret);

%'qFromx'
%printVec(shape(:,7))
%[iter Rval]
%[q dqdt]

%fid = fopen('qfx.txt','a');
%fprintf(fid,'%6d %+5.6e\n',iter,Rval);
%fclose(fid);