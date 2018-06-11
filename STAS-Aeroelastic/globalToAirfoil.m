function [wa,Dy,Tag,dTag] = globalToAirfoil (wg,Tas,TsB,TBB0,TB0g,dTsB,dTBB0)
%
% Computes w^a and the linear transform of a vector w^g.
%
%   States:           y vector:         u vector:
%                     qB
%                     qn1
%                     qn2
%                     wg
%
% Version:        Changes:
% --------        -------------
% 20.12.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 20.12.2017      Checked by complex step using wa output.
%
% Inputs:
% -------
% wg              : Vector, 3.
% Tas             : Transform from airfoil to section, constant.
% TsB,dTsB        : 3-by-3, 3-by-3*12, section-to-body and derivatives
%                   wrt qn1 and qn2.
% TBB0,dTBB0      : 3-by-3, 3-by-3*3, body CS derivatives wrt qB.
%
% Outputs:
% --------
% wa              : Vector, 3, airfoil coords.
% Dy              : State-space matrix.  dwa = Dy*[dqB,dqn1,dqn2,dwg].

Tsa = Tas.';
TBg = TB0g*TBB0;

Tag = TBg*TsB*Tas;
Tga = Tag.';
wa = Tga*wg;

Dy = zeros(3,21);
dTag = zeros(3,3*18);

Dy(:,19:21) = Tga;

for jj = 1:12
   jc3a = 3*(jj-1);
   jc3b = 3*(jj+6-1);
   dTag(:,jc3b+[1:3]) = TBg*dTsB(:,jc3a+[1:3])*Tas;
   Dy(:,jj+6) = (dTag(:,jc3b+[1:3]).')*wg;
end
for jj = 1:3  % Phi of ref node, indices 4:6 of Dy.
   jc3a = 3*(jj-1);
   jc3b = 3*(jj+3-1);
   dTag(:,jc3b+[1:3]) = TB0g*dTBB0(:,jc3a+[1:3])*TsB*Tas;
   Dy(:,jj+3) = (dTag(:,jc3b+[1:3]).')*wg;
end

%{
'---globalToAirfoil---'
wg
Tas
TsB
TBB0
TB0g
dTsB
dTBB0
'-----'
TBg
[Tag Tga]
wa
'---end gTA---'
%}
