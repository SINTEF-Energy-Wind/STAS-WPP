function dTsB = derivElementCS (qn1,qn2,Pn1,Pn2,TsB)
%
% This computes the derivatives of the T_s^B transform with
% respect to the input nodal positions and orientations. 
%
% Version:        Changes:
% --------        -------------
% 06.12.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 06.12.2017      Verified by complex step derivatives.
%
% Inputs:
% -------
% qn1, qn2        : Nodal DOF vectors.
% Pn1, Pn2        : Undeformed nodal vectors.
% TsB             : Transform, 3-by-3.
%
% Outputs:
% --------
% dTsB            : dTsB/dr and dTsB/dth.

dTsB = zeros(3,3*12);

Ts0Bjm1 = TFromTheta (Pn1(4:6));
Ts0Bj   = TFromTheta (Pn2(4:6));
[Tnn0k,dTnn0k] = dTdth (qn1(4:6));
[Tnn0k1,dTnn0k1] = dTdth (qn2(4:6));

Tn0ks0j = (Ts0Bj.')*Ts0Bjm1;
Teks0j  = Tn0ks0j*Tnn0k*(Tn0ks0j.');
Tek1s0j = Tnn0k1;

xs = TsB(:,1);
ys = TsB(:,2);
zs = TsB(:,3);

ysp = 0.5*Ts0Bj*(Teks0j(:,2) + Tek1s0j(:,2));

xn1 = zeros(6,1);
xn2 = zeros(6,1);
xn1(1:3) = Pn1(1:3) + qn1(1:3);
xn1(4:6) = qn1(4:6);
xn2(1:3) = Pn2(1:3) + qn2(1:3);
xn2(4:6) = qn2(4:6);

dxdr  = -dambnorm (xn1(1:3),xn2(1:3));
dydx  = -daxbnorm (xs,zs);
dydz  =  daxbnorm (zs,xs);
dzdx  =  daxbnorm (xs,ysp);
dzdyp = -daxbnorm (ysp,xs);

% Derivatives wrt xn1 r's.
for jj = 1:3
   jc3 = 3*(jj-1);
   dTsB(:,jc3+[1:3]) = [dxdr(:,jj)                             ...
                       dydx*dxdr(:,jj) + dydz*dzdx*dxdr(:,jj) ...
                       dzdx*dxdr(:,jj)];
end

% Derivatives wrt xn1 th's.
for jj = 1:3

   j3  = 3*(jj-1);
   jc3 = 3*(jj+3-1);

   dyspdth = 0.5*Ts0Bj*Tn0ks0j*dTnn0k(:,j3+[1:3])*(Tn0ks0j(2,:).');

   vec = dzdyp*dyspdth;
   dTsB(:,jc3+[1:3]) = [zeros(3,1) dydz*vec vec];

end

% Derivatives wrt xn2 r's.
dxdr = -dxdr;

for jj = 1:3
   jc3 = 3*(jj+6-1);
   dTsB(:,jc3+[1:3]) = [dxdr(:,jj)                             ...
                        dydx*dxdr(:,jj) + dydz*dzdx*dxdr(:,jj) ...
                        dzdx*dxdr(:,jj)];
end

% Derivatives wrt xn2 th's.
for jj = 1:3

   j3 = 3*(jj-1);
   jc3 = 3*(jj+9-1);

   dyspdth = 0.5*Ts0Bj*dTnn0k1(:,j3+2);

   vec = dzdyp*dyspdth;
   dTsB(:,jc3+[1:3]) = [zeros(3,1) dydz*vec vec];

end



