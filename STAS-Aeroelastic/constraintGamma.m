function [G,dLdq] = constraintGamma (q,dqdt,P,L,Lambda,Tbh, ... 
                                     ret,slv,idofs,idofm)
%
% The implementation of constraints in the equations of motion includes
% a term representing the centripetal constraint forces.  This term has
% the form Fc = Lambda^(-1)*M*Gamma, where Gamma includes a term
% -Ls^(-1)*dqk/dt*dL/dqk*dq/dt.
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
% q               : Nodal DOFs augmented with yaw, azimuth, front bearing
%                   axial offset, and three blade pitch DOFs.
% dqdt            : Rate of change of nodal DOFs.
% P               : Nodal offsets and initial body orientation.
% L               : Constraint equation matrix from constraints.m.
% Lambda          : Partitioned and inverted constraint matrix.
% Tbh             : Blade-to-hub coordinate transform.
% idofs,idofm     : Joint DOF references.
%
% Outputs:
% --------
% Gamma          : Partitioned and inverted constraint force matrix.
% dLdq           : Derivatives of perturbed constraint equations.

Neq = size(L,1);
Ndj = size(L,2);
Nret = size(ret,1);
Nslv = size(slv,1);

gdofs = getgdofs (idofs,idofm,Ndj);

% Derivative of constraints wrt nodal DOFs.
dLdq = derivConstraints (q,P,Tbh,idofs,idofm);

rsind = [ret;slv];

% Get the time derivative of L.
dLdt = spalloc (Neq,Ndj,Neq*Ndj);
for id = 1:size(gdofs,1)
   icN = Ndj*(id-1);
   dLdt = dLdt + dLdq(:,icN+rsind)*dqdt(gdofs(id));
end

Ls = L(:,slv);

G = sparse ([zeros(Nret,Nret);-(Ls\(dLdt*Lambda))]);



