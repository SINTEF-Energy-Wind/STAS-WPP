function [P,dPdq] = internalLoads (qn1,qn2,Pn1,Pn2,kes)
%
% In the postprocessing phase, compute internal loads and their TFs
% based on input nodal displacements and their TFs.
%
% Version:        Changes:
% --------        -------------
% 30.10.2019      Original code.
%
% Version:        Verification:
% --------        -------------
% 30.10.2019      Derivatives verified with complex step.
%
% Inputs:
% -------
% qn1, qn2        : Nodal DOF vectors.
% Pn1, Pn2        : Undeformed nodal vectors.
% kes             : Element stiffness matrix in section coordinates.
%
% Outputs:
% --------
% P               : Internal load vector, 12-by-1.
% dPdq            : Sensitivities wrt the input q's.


[xe,TsB] = elementCSFromNodes (qn1,qn2,Pn1,Pn2);
dTsB = derivElementCS (qn1,qn2,Pn1,Pn2,TsB);
mu = getMu (qn1,qn2,Pn1,Pn2,TsB);
dmu = dmudq (qn1,qn2,Pn1,Pn2,TsB,dTsB);

T = zeros(6,6);
for jj = 1:2
   iref = 3*(jj-1);
   T(iref+[1:3],iref+[1:3]) = TsB;
end

P = -T*kes(1:6,:)*mu;

dPdq = zeros(6,12);
for iq = 1:12
   ir3 = 3*(iq-1);
   dPdq(:,iq) = -T*kes(1:6,:)*dmu(:,iq) ...
      - [dTsB(:,ir3+[1:3]), zeros(3,3);zeros(3,3), dTsB(:,ir3+[1:3])]*kes(1:6,:)*mu;
end

%[qn1 qn2 Pn1 Pn2]
%[mu kes*mu [T, zeros(6,6);zeros(6,6), T]*kes*mu]



%{
% From buildK.m, MegaRoller FEM.

[xe,TsB] = elementCSFromNodes (qn1,qn2,Pn1,Pn2);
dTsB     = derivElementCS (qn1,qn2,Pn1,Pn2,TsB);
d2TsB    = secondDerivElementCS (qn1,qn2,Pn1,Pn2,TsB);
mu       = getMu (qn1,qn2,Pn1,Pn2,TsB);
dmu      = dmudq (qn1,qn2,Pn1,Pn2,TsB,dTsB);
d2mu     = d2mudq2 (qn1,qn2,Pn1,Pn2,TsB,dTsB,d2TsB);

kv = zeros(18,1);
for jj = 1:12
   kv(jj+6) = kv(jj+6) + (dmu(:,jj).')*kes*mu;
end

dkv = zeros(18,18);
mumu = mu*(mu.');
for ii = 1:12

   ic12 = 12*(ii-1);

   dmuimu = dmu(:,ii)*(mu.');

   for kk = 1:12

      % Verified.
      dkv(ii+6,kk+6) = dkv(ii+6,kk+6)                            ...
                     + (d2mu(:,ic12+kk).')*kes*mu ...
                     + (dmu(:,ii).')*kes*dmu(:,kk);

   end

end

P = kv(7:12);
dPdq = dkv(7:12,7:18);
%}