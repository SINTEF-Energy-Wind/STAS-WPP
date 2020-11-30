function [qB0,Pn0_B,Ts0_B,TB0_g] = ...
                    undeformedPosition (Pin,yaw,tilt,azimuth,cone,pitch,edx, ...
                                        idofs,idofm,inods,inodm)
%
% Computes the undeformed nodal positions and orientations.
%
% Version:        Changes:
% --------        -------------
% 02.10.2017      Original code.
% 27.11.2017      Updated to include joint angles at the end of qB0, Pn0_B.
%
% Version:        Verification:
% --------        -------------
% 02.10.2017      Visual check of outputs.
% 27.11.2017      Visual check of outputs.
%
% Inputs:
% -------
% Pin             : This is the vector of nodal coordinates that comes
%                   out of STASTurbine.  It contains pnB_B for each
%                   node, and then the structural twist angle.
% yaw...pitch     : yaw, azimuth, and pitch (1-3) are joint angles.
%                   tilt is shaft tilt, and cone is blade cone.
% edx             : Redundant DOF representing axial displacement of the
%                   driveshaft front bearing node with respect to the
%                   master node on the nacelle at this location.
% idofs...inodm   : DOF and node references for slave and master nodes.
%
% Outputs:
% --------
% qB0             : Since qB0 represents the deformation with respect
%                   to the undeformed shape, this is initialized to 
%                   zero.
% Pn0_B           : For body reference nodes: this gives OBg_g and
%                   the exponential map parameters Theta that describe
%                   the orientation of the body in the global CS.  For
%                   other nodes, this gives pnB_B and theta of the 
%                   associated element section CS.
% Ts0_B           : Transforms from element section to body coordinates.
% TB0_g           : Transforms from body to global coordinates.

Ndof = size(Pin,1);

qB0   = zeros(Ndof+6,1);
Pn0_B = zeros(Ndof+6,1);

% Initialize Pn0_B with Pin.
jj = [[1:6:Ndof-5] [2:6:Ndof-4] [3:6:Ndof-3]].';
Pn0_B(jj) = Pin(jj);

% Joint angles.
Pn0_B(Ndof+1) = yaw;
Pn0_B(Ndof+2) = azimuth;
Pn0_B(Ndof+3) = edx;
Pn0_B(Ndof+4) = pitch(1);
Pn0_B(Ndof+5) = pitch(2);
Pn0_B(Ndof+6) = pitch(3);

[Tn_y,Th_d,Tb_h] = basicTransforms (tilt,cone);
TB0_g = undeformedBodyToGlobal (yaw,tilt,azimuth,cone,pitch);

% Foundation and tower.
Pn0_B(idofs(1)+[1:3]) = zeros(3,1);

T = [0  1  0; ...
     0  0  1; ...
     1  0  0];
Pn0_B(idofs(1)+[4:6]) = thetaFromT (TB0_g(:,1:3));

vT  = thetaFromT (T);

for inod = inods(1)+1:inodm(2)
   ir3 = 3*(inod-1);
   ir6 = 6*(inod-1);
   Pn0_B(ir6+[4:6]) = vT;
end

Pn0_B(idofs(2)+[1:3]) = Pn0_B(idofs(1)+[1:3]) ...
                      + TB0_g(:,1:3)*Pn0_B(idofm(2)+[1:3]);
Pn0_B(idofs(2)+[4:6]) = thetaFromT (TB0_g(:,4:6));

for inod = inods(2)+1:inodm(3)
   ir3 = 3*(inod-1);
   ir6 = 6*(inod-1);
   Pn0_B(ir6+[4:6]) = vT;
end

% Nacelle.
Pn0_B(idofs(3)+[1:3]) = Pn0_B(idofs(2)+[1:3]) ...
                      + TB0_g(:,4:6)*Pn0_B(idofm(3)+[1:3]);
Pn0_B(idofs(3)+[4:6]) = thetaFromT (TB0_g(:,7:9));

nod = [inods(3)+1:inods(4)].';
for inod = 1:size(nod,1)-1
   ir6 = idofs(3) + 6*(inod-1);
   if (inod == 1)
      xs = Pn0_B(ir6+[7:9]);
   else
      xs = Pn0_B(ir6+[7:9]) - Pn0_B(ir6+[1:3]);
   end
   xs = xs/sqrt(xs.'*xs);
   ys = [0;-1;0];
   zs = cross(xs,ys);
   zs = zs/sqrt(zs.'*zs);
   T  = [xs ys zs];
   Pn0_B(ir6+[10:12]) = thetaFromT (T);
end

% Driveshaft.
T = [0 -1  0; ...
     0  0  1; ...
    -1  0  0];
Pn0_B(idofs(4)+[1:3]) = Pn0_B(idofs(3)+[1:3]) ...
                      + TB0_g(:,7:9)*Pn0_B(idofm(4)+[1:3]);
Pn0_B(idofs(4)+[4:6]) = thetaFromT (TB0_g(:,10:12));

for inod = inods(4)+1:inodm(6)-1
   ir6 = 6*(inod-1);
   Pn0_B(ir6+[4:6]) = thetaFromT (T);
end
for jj = 1:3
   inod = inodm(5+jj);
   ir3 = 3*(jj-1);
   ir6 = 6*(inod-1);
   Pn0_B(ir6+[4:6]) = thetaFromT (Th_d(:,ir3+[1:3]));
end

% Blades.
Nnb = inods(7) - inods(6);
for ib = 1:3
   Pn0_B(idofs(5+ib)+[1:3]) = Pn0_B(idofs(4)+[1:3]) ...
                            + TB0_g(:,10:12)*Pn0_B(idofm(5+ib)+[1:3]);
   Pn0_B(idofs(5+ib)+[4:6]) = thetaFromT (TB0_g(:,12+3*(ib-1)+[1:3]));

   [TB0_s,Ts0_B] = ...
         bladeSectionTransformsFromNodes (Pin(idofs(5+ib)+[1:6*Nnb]));

   for jj = 1:Nnb-1
      inod = inods(5+ib) + jj;
      ir3 = 3*(jj-1);
      ir6 = 6*(inod-1);
      Pn0_B(ir6+[4:6]) = thetaFromT (Ts0_B(:,ir3+[1:3]));
   end

end
