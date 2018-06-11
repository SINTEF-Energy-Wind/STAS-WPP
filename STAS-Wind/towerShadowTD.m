function Vgt = towerShadowTD (xeg,xhg,xyg,Vg,Tyg,Dt)
%
% Computes the upwind tower shadow effect for a set of input blade
% elements.  Inputs are in an arbitrary "global" coordinate system,
% with appropriate transforms.  This doesn't necessarily need to be
% the true global CS, it could be the rotorplane CS, for instance.
%
% The wind coordinate system has the x axis pointing downwind, the
% z axis along the centerline of the tower, and the y axis orthogonal.
% The tower shadow effect is computed by 2D potential flow theory
% using an equivalent position of the element relative to a 2D "slice"
% of the tower at the appropriate elevation.
%
% Version:        Changes:
% --------        -------------
% 04.03.2018      Original code.
%
% Version:        Verification:
% --------        -------------
% 04.03.2018      Checked some example cases.
%
% Inputs:
% -------
% xeg             : Element position in global coordinates.
% xhg             : Hub (rotorplane center) coordinates.
% xyg             : Yaw bearing coordinates.
% Vg              : Incoming xyz windspeed in global coordinates.  
%                   Should include induction, namely based on Vrxyz
%                   in BEMNL.
% Tyg             : Transform from yaw to global coordinates.
%                   We'll assume that the local deformation of the
%                   tower does not appreciably affect tower shadow.
% Dt              : Tower diameter at height of each element's annulus
%                   (size Nel).
%
% Outputs:
% --------
% Vgt             : Vg modified to account for tower shadow.

Nel = size(xeg,1)/3;

Tgy = Tyg.';

i3a = [1:3:3*Nel-2].';
i3b = [2:3:3*Nel-1].';
i3c = [3:3:3*Nel].';

% Rotorplane and yaw coordinates of the elements, expressed in
% the global frame.
rg = zeros(3*Nel,1);
dg = zeros(3*Nel,1);
rg(i3a) = xeg(i3a) - xhg(1);
rg(i3b) = xeg(i3b) - xhg(2);
rg(i3c) = xeg(i3c) - xhg(3);
dg(i3a) = xeg(i3a) - xyg(1);
dg(i3b) = xeg(i3b) - xyg(2);
dg(i3c) = xeg(i3c) - xyg(3);

% Transform to the yaw frame.
ry = zeros(3*Nel,1);
dy = zeros(3*Nel,1);
Vy = zeros(3*Nel,1);
for iel = 1:Nel
   ic3 = 3*(iel-1);
   ry(ic3+[1:3]) = Tgy*rg(ic3+[1:3]);
   dy(ic3+[1:3]) = Tgy*dg(ic3+[1:3]);
   Vy(ic3+[1:3]) = Tgy*Vg(ic3+[1:3]);
end

% Towerline offset vector.
sy = zeros(3*Nel,1);
sy(i3a) = dy(i3a);
sy(i3b) = ry(i3b);
sy(i3c) = sqrt(ry(i3b).^2 + ry(i3c).^2) + ry(i3c);

% Effective 2D coordinates.
xhat = sy(i3a);
sgn = sign(sy(i3b));
sgn(abs(sgn)<1) = 1;  % Need +/- 1, never zero.
yhat = sgn.*sqrt(sy(i3b).^2 + sy(i3c).^2);
%yhat = sqrt(sy(i3b).^2 + sy(i3c).^2);

% Transform to wind coordinates.
wang = atan2c (Vy(i3b),Vy(i3a));
cw = cos(wang);
sw = sin(wang);
x =  xhat.*cw + yhat.*sw;
y = -xhat.*sw + yhat.*cw;

% Velocity perturbation factor.
x2 = x.^2;
y2 = y.^2;
fact = -((0.5*Dt).^2).*(x2 - y2)./((x2 + y2).^2);
a = 1 + fact;
cw2 = cw.*cw;
sw2 = sw.*sw;
csw = cw.*sw;

% Reduced velocity.
Vyt = zeros(3*Nel,1);
Vyt(i3a) = Vy(i3a).*(a.*cw2 + sw2) + Vy(i3b).*csw.*(a-1);
Vyt(i3b) = Vy(i3a).*csw.*(a-1) + Vy(i3b).*(a.*sw2 + cw2);
Vyt(i3c) = Vy(i3c);

Vgt = zeros(3*Nel,1);
for iel = 1:Nel
   ic3 = 3*(iel-1);
   Vgt(ic3+[1:3]) = Tyg*Vyt(ic3+[1:3]);
end

%{
pang = atan2c (-ry(i3b),-ry(i3c))/pi;
'rg'
[pang rg(i3a) rg(i3b) rg(i3c)]
'dg'
[pang dg(i3a) dg(i3b) dg(i3c)]
'Vy'
[pang Vy(i3a) Vy(i3b) Vy(i3c)]
'sy'
[pang sy(i3a) sy(i3b) sy(i3c)]
'xhat yhat x y'
[pang xhat yhat x y]
'fact'
[pang fact]
%}

%{
fid = fopen('shad.txt','w');
for ip = 1:Nel
   ic3 = 3*(ip-1);
%   fprintf(fid,'%+5.6e %+5.6e %+5.6e %+5.6e %+5.6e\n', ...
%           xhat(ip),yhat(ip),x(ip),y(ip),fact(ip));
   fprintf(fid,'%+5.6e %+5.6e %+5.6e %+5.6e %+5.6e %+5.6e %+5.6e %+5.6e\n', ...
           atan2(-ry(ic3+2),-ry(ic3+3))/pi,fact(ip), ...
           Vyt(ic3+1),Vyt(ic3+2),Vgt(ic3+1),Vgt(ic3+2),Vgt(ic3+3), ...
           sqrt(Vgt(ic3+1)^2+Vgt(ic3+2)^2+Vgt(ic3+3)^2));
end
fclose('all');
%}
