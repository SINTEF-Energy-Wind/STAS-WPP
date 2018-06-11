function [xe,TsB] = elementCSFromNodes (qn1,qn2,Pn1,Pn2)
%
% This computes the element section coordinate system origin and
% orientation from the values at the adjacent nodes.  Here qn1 holds
% the d and theta values associated with node k.  The d's are given
% in body coordinates.  The theta's specify the T_n^n0 transform.  At
% node k, the "n0" coordinate system is aligned with the section CS 
% of the previous element j-1.  (For element 1 this can be taken to
% be equal to element 1's section CS.) The qn2 input holds the values
% d and theta for node k+1, where the "n0" coordinate system is now
% that of element j.  Pn1 and Pn2 contain rotation parameters that
% specify the T_s0^B transforms for elements j-1 and j, respectively.
%
% Version:        Changes:
% --------        -------------
% 06.12.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 06.12.2017      Selected cases checked and documented in report.
%                 Also 45-degree bend case.
%
% Inputs:
% -------
% qn1, qn2        : Nodal vectors.  The position parts are dn^B,
%                   and the rotations code for T_n^n0.
% Pn1,Pn2         : 6*Nel vector describing undeformed nodal position
%                   and orientation of nodes.  By definition, Pn1 and 
%                   Pn2 contain the orientation of the undeformed
%                   section coordinate systems, giving T_s0^B.
%
% Outputs:
% --------
% xe              : Element position wrt body CS.
% TsB             : Transform from element section to body CS.

%'elementCSFromNodes'

Ts0Bjm1 = TFromTheta (Pn1(4:6));
Ts0Bj   = TFromTheta (Pn2(4:6));
Tnn0k   = TFromTheta (qn1(4:6));
Tnn0k1  = TFromTheta (qn2(4:6));

Tn0ks0j = (Ts0Bj.')*Ts0Bjm1;
Teks0j  = Tn0ks0j*Tnn0k*(Tn0ks0j.');
Tek1s0j = Tnn0k1;

xe = zeros(6,1);
xe(1:3) = 0.5*(qn1(1:3) + qn2(1:3) + Pn1(1:3) + Pn2(1:3));

dr = Pn2(1:3) + qn2(1:3) - Pn1(1:3) - qn1(1:3);
TsB(:,1) = dr/sqrt((dr.')*dr);

ysp = 0.5*Ts0Bj*(Teks0j(:,2) + Tek1s0j(:,2));

vec = cross(TsB(:,1),ysp);
TsB(:,3) = vec/sqrt((vec.')*vec);

vec = cross(TsB(:,3),TsB(:,1));
TsB(:,2) = vec/sqrt((vec.')*vec);

xe(4:6) = thetaFromT (TsB);


%{
'--------- TsB ----------'

[Pn1 qn1 Pn2 qn2]
Ts0Bjm1
Ts0Bj
Tn0ks0j
Teks0j
Tek1s0j
TsB

TsB-Ts0Bj

'------- end TsB --------'
%}

