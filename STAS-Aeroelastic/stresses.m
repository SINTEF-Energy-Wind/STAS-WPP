function [s,dsdq] = stresses (qn1,qn2,Pn1,Pn2,ry,rz,EE)
%
% In the postprocessing phase, compute stresses and their gradients
% based on input nodal displacements.
%
% Version:        Changes:
% --------        -------------
% 10.09.2020      Original code.
%
% Version:        Verification:
% --------        -------------
% 10.09.2020      
%
% Inputs:
% -------
% qn1, qn2        : Nodal DOF vectors.
% Pn1, Pn2        : Undeformed nodal vectors.
% ry, rz          : Coordinates of points about the section at
%                   which to compute the stress.                   
% EE              : 3-by-3 elasticity matrix in section coordinates,
%                   for each point.  Packed as 3-by-3*Nsp.
%
% Outputs:
% --------
% s               : Stresses.
% dsdq            : Sensitivities wrt the input q's.

Nsp = size(ry,1);

[xe,TsB] = elementCSFromNodes (qn1,qn2,Pn1,Pn2);
dTsB = derivElementCS (qn1,qn2,Pn1,Pn2,TsB);
mu = getMu (qn1,qn2,Pn1,Pn2,TsB);
dmu = dmudq (qn1,qn2,Pn1,Pn2,TsB,dTsB);

L = sqrt((Pn2(1)-Pn1(1))^2 + (Pn2(2)-Pn1(2))^2 + (Pn2(3)-Pn1(3))^2);
Li = 1/L;
Li2 = Li^2;

dSmat = [-Li,   0,    0,   0,   0,   0,   Li,  0,    0,   0,   0,   0; ...
          0,    0,    0,   0,   0,   1,   0,   0,    0,   0,   0,   0; ...
          0,    0,    0,   0,  -1,   0,   0,   0,    0,   0,   0,   0; ...
          0,    0,    0,  -Li,  0,   0,   0,   0,    0,   Li,  0,   0; ...
          0,    0,  6*Li2, 0, -4*Li, 0,   0,   0, -6*Li2, 0, -2*Li, 0; ...
          0, -6*Li2,  0,   0,   0, -4*Li, 0, 6*Li2,  0,   0,   0, -2*Li];

s = zeros(3,Nsp);
dsdq = zeros(3,12*Nsp);
for ip = 1:Nsp

   ic3 = 3*(ip-1);
   ic12 = 12*(ip-1);

   Dmat = [1, 0, 0,    0,    rz(ip), -ry(ip); ...
           0, 0, 0, -rz(ip),   0,       0;    ...
           0, 0, 0,  ry(ip),   0,       0];

   EDS = EE(:,ic3+[1:3])*Dmat*dSmat;
   s(:,ip) = EDS*mu;

   dsdq(:,ic12+[1:12]) = EDS*dmu;

end

