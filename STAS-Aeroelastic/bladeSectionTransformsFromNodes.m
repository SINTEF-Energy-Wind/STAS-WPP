function [TB_s,Ts_B] = bladeSectionTransformsFromNodes (Pin)
%
% This computes body section transforms from dx,y,z and twist
% inputs contained in the Pin vector.
%
% Version:        Changes:
% --------        -------------
% 17.09.2017      Original code
%
% Version:        Verification:
% --------        -------------
% 17.09.2017      Tested for sample x,y,z,twist inputs.
%
% Inputs:
% -------
% Pin             : Nodal positions and twist-bend-sweep rotation.
%                   Of the rotations, only twist is used.  Bend and
%                   sweep are defined by nodal positions.
%
% Outputs:
% --------
% TB_s,Ts_B

%'bladeSectionTransformsFromNodes'

Nnod = size(Pin,1)/6;
Nel  = Nnod - 1;

TB_s = zeros(3,3*Nel);
Ts_B = zeros(3,3*Nel);

dr  = zeros(3*Nel,1);
L   = zeros(Nel,1);
for idof = 1:3
   ii = [idof:3:3*(Nel-1)+idof];
   jj = [idof:6:6*(Nnod-2)+idof];
   kk = [idof+6:6:6*(Nnod-1)+idof];

   dr(ii) = Pin(kk) - Pin(jj);

   L(:,1) = L(:,1) + dr(ii).^2;

end
L = sqrt(L);

% Do not average.  The element values are stored in the outboard node.
%xie = 0.5*(Pin(4:6:6*Nel-2) + Pin(10:6:6*Nnod-2));
xie = Pin(10:6:6*Nnod-2);

c1 = cos(xie);              % Twist (Rotation about body x axis).
s1 = sin(xie);

th2 = asin(-dr(3:3:3*Nel)./L);
th3 = atan2c(dr(2:3:3*Nel-1),dr(1:3:3*Nel-2));

c2 = cos(th2);              % Rotation about body y axis.
s2 = sin(th2);
c3 = cos(th3);              % Rotation about body z axis.
s3 = sin(th3);

for iel = 1:Nel

   ic = 3*(iel-1);

   TB_s(:,ic+[1:3]) = [1 0 0;0 c1(iel) s1(iel);0 -s1(iel) c1(iel)] ...
                    * [c2(iel) 0 -s2(iel);0 1 0;s2(iel) 0 c2(iel)] ...
                    * [c3(iel) s3(iel) 0;-s3(iel) c3(iel) 0;0 0 1];

   Ts_B(:,ic+[1:3]) = TB_s(:,ic+[1:3]).';

%iel
%Ts_B(:,ic+[1:3])
%dr = Pin(6*(iel-1)+[7:9]) - Pin(6*(iel-1)+[1:3]);
%[Pin(6*(iel-1)+[1:3]) Pin(6*(iel-1)+[7:9]) dr]
%dr/sqrt(dr.'*dr)

%{
'BladeSectionTransformsFromNodes'
iel
[c3(iel) s3(iel) 0;-s3(iel) c3(iel) 0;0 0 1]
[c2(iel) 0 -s2(iel);0 1 0;s2(iel) 0 c2(iel)]
[1 0 0;0 c1(iel) s1(iel);0 -s1(iel) c1(iel)]
Ts_B(:,ic+[1:3])
%}

end


