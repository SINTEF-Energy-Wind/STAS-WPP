function [blxdof,bludof,blydof] = ...
             MBCindices (Nxa,Nxp,Neta,Ndj,Neb,imdofs,idofs)

Nxab = Nxa/3;
Nxpb = Nxp/3;

Nxs  = 2*Neta;
Nus  = Ndj + Neta;

Nx1  = Nxs + Nxa;

[b1eta,b2eta,b3eta] = MBCindices_Nmd (Neta,imdofs);
[b1q,b2q,b3q]       = MBCindices_Ndj (Ndj,idofs);

Nblxdof = 2*size(b1eta,1) + Nxab + Nxpb;
blxdof = zeros(Nblxdof,3);
blxdof(:,1) = [b1eta.', Neta+b1eta.', Nxs+[1:Nxab],        Nx1+[1:Nxpb]].';
blxdof(:,2) = [b2eta.', Neta+b2eta.', Nxs+Nxab+[1:Nxab],   Nx1+Nxpb+[1:Nxpb]].';
blxdof(:,3) = [b3eta.', Neta+b3eta.', Nxs+2*Nxab+[1:Nxab], Nx1+2*Nxpb+[1:Nxpb]].';

Nbludof = size(b1q,1) + size(b1eta,1) + 3*Neb;
bludof = zeros(Nbludof,3);
bludof(:,1) = [b1q.', Ndj+b1eta.', Nus+[1:3*Neb]].';
bludof(:,2) = [b2q.', Ndj+b2eta.', Nus+3*Neb+[1:3*Neb]].';
bludof(:,3) = [b3q.', Ndj+b3eta.', Nus+6*Neb+[1:3*Neb]].';

Nblydof = 4*size(b1q,1) + 1 + 1;
blydof = zeros(Nblydof,3);
blydof(:,1) = [b1q.', Ndj+b1q.', 2*Ndj+b1q.', 3*Ndj+b1q.', 4*Ndj+1, 4*Ndj+7+1].';
blydof(:,2) = [b2q.', Ndj+b2q.', 2*Ndj+b2q.', 3*Ndj+b2q.', 4*Ndj+2, 4*Ndj+7+2].';
blydof(:,3) = [b3q.', Ndj+b3q.', 2*Ndj+b3q.', 3*Ndj+b3q.', 4*Ndj+3, 4*Ndj+7+3].';