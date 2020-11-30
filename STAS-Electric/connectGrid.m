function [dxdt,A,B] = connectGrid (linFlag,x,u,p,nod)
%
% State-space model of a generic AC grid, cables connecting buses.
%
%             States:           u vector:
%   Overall:
%                               we          (Reference synchronous gen.) 
%   For each cable segment:
%             iLd,q      
%
%   For each bus:
%             vd,q
%
%   For each bus:
%             ird,q             ied,q       (External current entering bus.)
%
% Version:        Changes:
% --------        -------------
% 08.03.2019      Original code.
%
% Version:        Verification:
% --------        -------------
% 08.03.2019      Derivatives verified by complex step.
%
% Inputs:
% -------
% p               : Ncab:  R     (Ohms)  Phase resistance in line j
%                   Ncab:  L     (H)     Phase inductance in line j
%                   Nnod:  C     (F)     Phase capacitance to ground at node k
%                   Nnod:  Rr    (Ohms)  Phase resistance to ground at node k.
%                   Nnod:  Lr    (H)     Phase inductance to ground at node k.
% nod             : Connectivity, 2-by-Nnod.
%
% Outputs:
% --------
%

Ncab = size(nod,2);
Nnod = max(max(nod));

Nx = 2*Ncab + 4*Nnod;
Nu = 1 + 2*Nnod;

ixi = 0;
ixv = 2*Ncab;
ixr = ixv + 2*Nnod;
iue = 1;

we = u(1);
R  = p(1:Ncab);
L  = p(Ncab+[1:Ncab]);
C  = p(2*Ncab+[1:Nnod]);
Rr = p(2*Ncab+Nnod+[1:Nnod]);
Lr = p(2*Ncab+2*Nnod+[1:Nnod]);

dxdt = zeros(Nx,1);

spn = [0 -1;1 0];

for icab = 1:Ncab

   ic2 = 2*(icab-1);

   if (L(icab) > 0)
      dxdt(ixi+ic2+[1:2]) = -(we*spn + (R(icab)/L(icab))*speye(2))*x(ixi+ic2+[1:2]) ...
                          + (x(ixv+2*nod(1,icab)+[-1:0])                            ...
                          -  x(ixv+2*nod(2,icab)+[-1:0]))/L(icab);
   end

end

i2a = ixi + [1:2:2*Ncab-1].';
i2b = ixi + [2:2:2*Ncab].';

for inod = 1:Nnod

   ic2 = 2*(inod-1);

   ind1 = (nod(1,:) == inod).';
   ind2 = (nod(2,:) == inod).';

   if (C(inod) > 0)
      dxdt(ixv+ic2+[1:2]) = -we*spn*x(ixv+ic2+[1:2])                ...
                          + (-[sum(x(i2a(ind1)));sum(x(i2b(ind1)))]  ...
                          +   [sum(x(i2a(ind2)));sum(x(i2b(ind2)))]  ...
                          +  u(iue+ic2+[1:2]))/C(inod);
   end

   if (Lr(inod) > 0)
      dxdt(ixr+ic2+[1:2]) = -(we*spn + (Rr(inod)/Lr(inod))*speye(2))*x(ixr+ic2+[1:2]) ...
                          + x(ixv+ic2+[1:2])/Lr(inod);
   end

end



if (linFlag == 1)

   A = spalloc(Nx,Nx,0.1*Nx*Nx);
   B = spalloc(Nx,Nu,0.1*Nx*Nu);

   for icab = 1:Ncab

      ic2 = 2*(icab-1);

      if (L(icab) > 0)
         ir = ixi + ic2 + [1:2];
         ic = ir;
         A(ir,ic) = -(we*spn + (R(icab)/L(icab))*speye(2));
         ic = ixv+2*nod(1,icab)+[-1:0];
         A(ir,ic) = A(ir,ic) + speye(2)/L(icab);
         ic = ixv+2*nod(2,icab)+[-1:0];
         A(ir,ic) = A(ir,ic) - speye(2)/L(icab);
         ic = 1;
         B(ir,ic) = -spn*x(ixi+ic2+[1:2]);
      end

   end

   for inod = 1:Nnod

      ic2 = 2*(inod-1);

      ind1 = (nod(1,:) == inod).';
      ind2 = (nod(2,:) == inod).';

      if (C(inod) > 0)
         ir = ixv+ic2+[1:2];
         ic = ir;
         A(ir,ic) = A(ir,ic) - we*spn;
         ir = ixv+ic2+1;
         ic = i2a(ind1);
         A(ir,ic) = A(ir,ic) - ones(1,sum(ind1))/C(inod);
         ic = i2a(ind2);
         A(ir,ic) = A(ir,ic) + ones(1,sum(ind2))/C(inod);
         ir = ixv+ic2+2;
         ic = i2b(ind1);
         A(ir,ic) = A(ir,ic) - ones(1,sum(ind1))/C(inod);
         ic = i2b(ind2);
         A(ir,ic) = A(ir,ic) + ones(1,sum(ind2))/C(inod);
         ir = ixv+ic2+[1:2];
         ic = iue+ic2+[1:2];
         B(ir,ic) = B(ir,ic) + speye(2)/C(inod);
         ic = 1;
         B(ir,ic) = B(ir,ic) - spn*x(ixv+ic2+[1:2]);

      end

      if (Lr(inod) > 0)
         ir = ixr+ic2+[1:2];
         ic = ir;
         A(ir,ic) = A(ir,ic) - we*spn - (Rr(inod)/Lr(inod))*speye(2);
         ic = ixv+ic2+[1:2];
         A(ir,ic) = A(ir,ic) + speye(2)/Lr(inod);
         ic = 1;
         B(ir,ic) = B(ir,ic) - spn*x(ixr+ic2+[1:2]);
      end

   end

else

   A = sparse(Nx,Nx);
   B = sparse(Nx,Nu);

end

