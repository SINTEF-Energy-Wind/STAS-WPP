function [x0,dq0dt,ret,slv] = initialDOFVelocity (s,q,P,Omega)
% 
% Compute the initial DOFs' rate of change based on a given deflected
% or undeflected shape and rotor speed.
%
% Version:        Changes:
% --------        -------------
% 19.03.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 19.03.2017      Checked for a simple sample case.
%
% Inputs:
% -------
% 
%
% Outputs:
% --------
% 

[idofs,idofm,inods,inodm,Ndof] = getDOFRefs (s);
Ndj = size(q,1);

% Retained and slave structural DOFs.
slv = slaveDOFs (idofs);
vec = [1:Ndj].';
[jnk,ret,jnk2] = partitionMatrix (vec,slv,[]);
Nret = size(ret,1);
Nslv = size(slv,1);

x0   = zeros(2*Nret,1);
x0(1:Nret) = q(ret);
x0(2*Nret-6+2) = Omega;
qg   = zeros(Nslv,1);
[q,dq0dt] = qFromx (x0,s,P,speye(Nret),ret,slv,qg);
