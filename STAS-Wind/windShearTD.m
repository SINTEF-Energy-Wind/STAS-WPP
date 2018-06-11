function Vgs = windShearTD (typ,xeg,xhg,Vg,param)
%
% Computes the wind shear effect.
%
% Version:        Changes:
% --------        -------------
% 05.03.2018      
%
% Version:        Verification:
% --------        -------------
% 05.03.2018      
%
% Inputs:
% -------
% typ             : = 1 for logarithmic, = 2 for power law.
% xeg             : Element coordinates.
% xhg             : Hub center coordinates or an alternate reference
%                   height.
% Vgh             : Nominal windspeed at the reference height.
% param           : Surface roughness in the event that typ = 1.
%                   Power exponent in the event that typ = 2.
%
% Outputs:
% --------
% Vgs             : Windspeed with shear.

Nel = size(xeg,1)/3;

i3a = [1:3:3*Nel-2].';
i3b = [2:3:3*Nel-1].';
i3c = [3:3:3*Nel].';

z = xeg(i3c);

H0 = xhg(3);

if (typ == 1)
   f = log(z/param)/log(H0/param);
elseif (typ == 2)
   f = (z/H0).^param;
end

Vgs = zeros(3*Nel,1);
Vgs(i3a) = Vg(i3a).*f;
Vgs(i3b) = Vg(i3b).*f;
Vgs(i3c) = Vg(i3c).*f; 

%{
fid = fopen('shear.txt','w');
for ip = 1:Nel
   ic3 = 3*(ip-1);
   fprintf(fid,'%+5.6e %+5.6e %+5.6e %+5.6e %+5.6e\n', ...
           z(ip),f(ip),Vgs(ic3+1),Vgs(ic3+2),Vgs(ic3+3));
end
fclose('all');
%}