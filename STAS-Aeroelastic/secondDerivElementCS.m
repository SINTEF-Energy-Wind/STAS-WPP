function d2TsB = secondDerivElementCS (qn1,qn2,Pn1,Pn2,TsB)
%
% Computes d^2T_s^B/dr^2 /drdth and /dth^2.
%
% The transform is a function of the positions and rotations of
% nodes 1 and 2.
%
% This initial implementation is intentionally inefficient, to aid
% debugging and verification.  A future implementation should take
% advantage of symmetries.
%
% Version:        Changes:
% --------        -------------
% 07.11.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 07.11.2017      Verified using complex step derivatives.
%
% Inputs:
% -------
% qn1, qn2        : Nodal DOF vectors.
% Pn1, Pn2        : Undeformed nodal vectors.
% TsB             : Transform, 3-by-3.
%
% Outputs:
% --------        (x3)
% d2TsB           : 1: d2TsB/(dr1x  dr1x)
%                   2: d2TsB/(dr1x  dr1y)
%                   ...
%                  12: d2TsB/(dr1x  dth2z)
%                  13: d2TsB/(dr1y  dr1x)
%                  14: d2TsB/(dr1y  dr1y)
%                   ...
%                  24: d2TsB/(dr1y  dth2z)
%                   ...
%                  61: d2TsB/(dth1z dr1x)
%                  62: d2TsB/(dth1z dr1y)
%                   ...
%                  72: d2TsB/(dth1z dth2z)
%                  73: d2TsB/(dr2x  dr1x)
%                   ...
%                 144: d2TsB/(dth2z dth2z)

d2TsB = zeros (3,3*12*12);

Ts0Bjm1 = TFromTheta (Pn1(4:6));
Ts0Bj   = TFromTheta (Pn2(4:6));
[Tnn0k,dTnn0k] = dTdth (qn1(4:6));
d2Tnn0k = d2Tdth2 (qn1(4:6),Tnn0k,dTnn0k);
[Tnn0k1,dTnn0k1] = dTdth (qn2(4:6));
d2Tnn0k1 = d2Tdth2 (qn2(4:6),Tnn0k1,dTnn0k1);

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

dxdr     =  dambnorm (xn1(1:3),xn2(1:3));  % CHECK SIGN.  This matches finite
                                           % difference, which means I must have
                                           % missed a sign elsewhere... but the
                                           % final answer is correct.
dydx     = -daxbnorm (xs,zs);
dydz     =  daxbnorm (zs,xs);
dzdx     =  daxbnorm (xs,ysp);
dzdyp    = -daxbnorm (ysp,xs);

d2xdr2   =  d2ambnorm (xn2(1:3),xn1(1:3));  % d^2x/dr2^2.  d^2x/(dr1 dr2) = -d^2x/dr2^2.
                                            % d^2x/dr1^2 = d^2x/dr2^2.  Verified.
d2ydx2   = -d2axbnorm (xs,zs);              % Verified.
d2ydz2   =  d2axbnorm (zs,xs);              % Verified.
d2zdx2   =  d2axbnorm (xs,ysp);             % Verified.
   
d2ydxdz  = -d2axbdab (xs,zs);               % Verified.
d2zdxdyp =  d2axbdab (xs,ysp);              % 

d2zdyp2  = -d2axbnorm (ysp,xs);             % Verified.

dzdth = zeros(3,6);
dypdth = zeros(3,6);                        % Verified.
for jj = 1:3
   j3 = 3*(jj-1);
   dypdth(:,jj)   = 0.5*Ts0Bj*Tn0ks0j*dTnn0k(:,j3+[1:3])*(Tn0ks0j(2,:).');
   dypdth(:,jj+3) = 0.5*Ts0Bj*dTnn0k1(:,j3+2);
   dzdth(:,jj)    = dzdyp*dypdth(:,jj);
   dzdth(:,jj+3)  = dzdyp*dypdth(:,jj+3);
end

% Cross-derivatives between qn1 and qn2 are zero.  Cross-derivatives of
% components within qn1 or qn2 may be non-zero.  [Verified.]
d2ypdth2 = zeros(3,6*6);
for ii = 1:3
   i9 = 9*(ii-1);
   ic6 = 6*(ii-1);
   for jj = 1:3
      j3 = 3*(jj-1);
      d2ypdth2(:,ic6+jj) = 0.5*Ts0Bj*Tn0ks0j*d2Tnn0k(:,i9+j3+[1:3])*(Tn0ks0j(2,:).');
   end
end
for ii = 1:3
   i9 = 9*(ii-1);
   ic6 = 6*(ii+3-1);
   for jj = 1:3
      j3 = 3*(jj-1);
      d2ypdth2(:,ic6+3+jj) = 0.5*Ts0Bj*d2Tnn0k1(:,i9+j3+2);
   end
end

% d2TsB/(dr1 dr1).
% The sign on dx/dr_i is negative.  The sign on dx/dr_j is negative.
% The sign on d2xdr2 is positive.
for ii = 1:3

   i36 = 36*(ii-1);  % r1
   i3  =  3*(ii-1);

   dxdri = -dxdr(:,ii);
   dzdri = dzdx*dxdri;

   for jj = 1:3

      j3 = 3*(jj-1);  % r1
      indx = i36 + j3;

      dxdrj = -dxdr(:,jj);
      d2xij = d2xdr2(:,i3+jj);
      dzdrj = dzdx*dxdrj;

      d2TsB(:,indx+1) = d2xij;

      d2TsB(:,indx+2) = dydx*d2xij;
      d2TsB(:,indx+2) = d2TsB(:,indx+2) + dydz*dzdx*d2xij;
      d2TsB(:,indx+3) = dzdx*d2xij;
      for pp = 1:3
         p3 = 3*(pp-1);
         dxdx  = dxdrj*dxdri(pp);
         d2zdxdx = d2zdx2(:,p3+[1:3])*dxdx;
         term1 = d2ydx2(:,p3+[1:3])*dxdx;
         term2 = d2ydz2(:,p3+[1:3])*dzdrj*dzdri(pp);
         term3 = dydz*d2zdxdx;
         term4 = d2ydxdz(:,p3+[1:3])*dzdrj*dxdri(pp);
         term5 = d2ydxdz(:,pp+[0 3 6])*dxdrj*dzdri(pp);
         d2TsB(:,indx+2) = d2TsB(:,indx+2) ...
                         + term1 + term2 + term3 + term4 + term5;
         d2TsB(:,indx+3) = d2TsB(:,indx+3) + d2zdxdx;
      end

   end      

end

% d2TsB/(dr1 dr2).
% The sign on dx/dr_i is negative.  The sign on dx/dr_j is positive.
% The sign on d2xdr2 is negative.
for ii = 1:3

   i36 = 36*(ii-1);  % r1
   i3  =  3*(ii-1);

   dxdri = -dxdr(:,ii);
   dzdri = dzdx*dxdri;

   for jj = 1:3

      j3 = 3*(jj+6-1);  % r2
      indx = i36 + j3;

      dxdrj = dxdr(:,jj);
      d2xij = -d2xdr2(:,i3+jj);
      dzdrj = dzdx*dxdrj;

      d2TsB(:,indx+1) = d2xij;

      d2TsB(:,indx+2) = dydx*d2xij;
      d2TsB(:,indx+2) = d2TsB(:,indx+2) + dydz*dzdx*d2xij;
      d2TsB(:,indx+3) = dzdx*d2xij;
      for pp = 1:3
         p3 = 3*(pp-1);
         dxdx  = dxdrj*dxdri(pp);
         d2zdxdx = d2zdx2(:,p3+[1:3])*dxdx;
         term1 = d2ydx2(:,p3+[1:3])*dxdx;
         term2 = d2ydz2(:,p3+[1:3])*dzdrj*dzdri(pp);
         term3 = dydz*d2zdxdx;
         term4 = d2ydxdz(:,p3+[1:3])*dzdrj*dxdri(pp);
         term5 = d2ydxdz(:,pp+[0 3 6])*dxdrj*dzdri(pp);
         d2TsB(:,indx+2) = d2TsB(:,indx+2) ...
                         + term1 + term2 + term3 + term4 + term5;
         d2TsB(:,indx+3) = d2TsB(:,indx+3) + d2zdxdx;
      end

   end      

end

% d2TsB/(dr2 dr1).
% The sign on dx/dr_i is positive.  The sign on dx/dr_j is negative.
% The sign on d2xdr2 is negative.
for ii = 1:3

   i36 = 36*(ii+6-1);  % r2
   i3  =  3*(ii-1);

   dxdri = dxdr(:,ii);
   dzdri = dzdx*dxdri;

   for jj = 1:3

      j3 = 3*(jj-1);  % r1
      indx = i36 + j3;

      dxdrj = -dxdr(:,jj);
      d2xij = -d2xdr2(:,i3+jj);
      dzdrj = dzdx*dxdrj;

      d2TsB(:,indx+1) = d2xij;

      d2TsB(:,indx+2) = dydx*d2xij;
      d2TsB(:,indx+2) = d2TsB(:,indx+2) + dydz*dzdx*d2xij;
      d2TsB(:,indx+3) = dzdx*d2xij;
      for pp = 1:3
         p3 = 3*(pp-1);
         dxdx  = dxdrj*dxdri(pp);
         d2zdxdx = d2zdx2(:,p3+[1:3])*dxdx;
         term1 = d2ydx2(:,p3+[1:3])*dxdx;
         term2 = d2ydz2(:,p3+[1:3])*dzdrj*dzdri(pp);
         term3 = dydz*d2zdxdx;
         term4 = d2ydxdz(:,p3+[1:3])*dzdrj*dxdri(pp);
         term5 = d2ydxdz(:,pp+[0 3 6])*dxdrj*dzdri(pp);
         d2TsB(:,indx+2) = d2TsB(:,indx+2) ...
                         + term1 + term2 + term3 + term4 + term5;
         d2TsB(:,indx+3) = d2TsB(:,indx+3) + d2zdxdx;
      end

   end      

end

% d2TsB/(dr2 dr2).
% The sign on dx/dr_i is positive.  The sign on dx/dr_j is positive.
% The sign on d2xdr2 is positive.
for ii = 1:3

   i36 = 36*(ii+6-1);  % r2
   i3  =  3*(ii-1);

   dxdri = dxdr(:,ii);
   dzdri = dzdx*dxdri;

   for jj = 1:3

      j3 = 3*(jj+6-1);  % r2
      indx = i36 + j3;

      dxdrj = dxdr(:,jj);
      d2xij = d2xdr2(:,i3+jj);
      dzdrj = dzdx*dxdrj;

      d2TsB(:,indx+1) = d2xij;

      d2TsB(:,indx+2) = dydx*d2xij;
      d2TsB(:,indx+2) = d2TsB(:,indx+2) + dydz*dzdx*d2xij;
      d2TsB(:,indx+3) = dzdx*d2xij;
      for pp = 1:3
         p3 = 3*(pp-1);
         dxdx  = dxdrj*dxdri(pp);
         d2zdxdx = d2zdx2(:,p3+[1:3])*dxdx;
         term1 = d2ydx2(:,p3+[1:3])*dxdx;
         term2 = d2ydz2(:,p3+[1:3])*dzdrj*dzdri(pp);
         term3 = dydz*d2zdxdx;
         term4 = d2ydxdz(:,p3+[1:3])*dzdrj*dxdri(pp);
         term5 = d2ydxdz(:,pp+[0 3 6])*dxdrj*dzdri(pp);
         d2TsB(:,indx+2) = d2TsB(:,indx+2) ...
                         + term1 + term2 + term3 + term4 + term5;
         d2TsB(:,indx+3) = d2TsB(:,indx+3) + d2zdxdx;
      end

   end      

end

% d2TsB/(dr1 dth1).
% The sign on dx/dr_i is negative.  
for ii = 1:3

   i36 = 36*(ii-1);  % r1

   dxdri = -dxdr(:,ii);
   dzdri = dzdx*dxdri;

   for jj = 1:3

      j3 = 3*(jj+3-1);  % th1
      indx = i36 + j3;

      for pp = 1:3
         p3 = 3*(pp-1);
         term1 = d2ydxdz(:,p3+[1:3])*dzdth(:,jj)*dxdri(pp);
         term2 = d2ydz2(:,p3+[1:3])*dzdth(:,jj)*dzdri(pp);
         term3 = dydz*d2zdxdyp(:,p3+[1:3])*dypdth(:,jj)*dxdri(pp);
         d2TsB(:,indx+2) = d2TsB(:,indx+2) + term1 + term2 + term3;
         d2TsB(:,indx+3) = d2TsB(:,indx+3) ...
                         + d2zdxdyp(:,p3+[1:3])*dypdth(:,jj)*dxdri(pp);
      end

   end      

end

% d2TsB/(dr1 dth2).
% The sign on dx/dr_i is negative.  
for ii = 1:3

   i36 = 36*(ii-1);  % r1

   dxdri = -dxdr(:,ii);
   dzdri = dzdx*dxdri;

   for jj = 1:3

      j3 = 3*(jj+9-1);  % th2
      indx = i36 + j3;

      for pp = 1:3
         p3 = 3*(pp-1);
         term1 = d2ydxdz(:,p3+[1:3])*dzdth(:,jj+3)*dxdri(pp);
         term2 = d2ydz2(:,p3+[1:3])*dzdth(:,jj+3)*dzdri(pp);
         term3 = dydz*d2zdxdyp(:,p3+[1:3])*dypdth(:,jj+3)*dxdri(pp);
         d2TsB(:,indx+2) = d2TsB(:,indx+2) + term1 + term2 + term3;
         d2TsB(:,indx+3) = d2TsB(:,indx+3) ...
                         + d2zdxdyp(:,p3+[1:3])*dypdth(:,jj+3)*dxdri(pp);
      end

   end      

end

% d2TsB/(dr2 dth1).
% The sign on dx/dr_i is positive.  
for ii = 1:3

   i36 = 36*(ii+6-1);  % r2

   dxdri = dxdr(:,ii);
   dzdri = dzdx*dxdri;

   for jj = 1:3

      j3 = 3*(jj+3-1);  % th1
      indx = i36 + j3;

      for pp = 1:3
         p3 = 3*(pp-1);
         term1 = d2ydxdz(:,p3+[1:3])*dzdth(:,jj)*dxdri(pp);
         term2 = d2ydz2(:,p3+[1:3])*dzdth(:,jj)*dzdri(pp);
         term3 = dydz*d2zdxdyp(:,p3+[1:3])*dypdth(:,jj)*dxdri(pp);
         d2TsB(:,indx+2) = d2TsB(:,indx+2) + term1 + term2 + term3;
         d2TsB(:,indx+3) = d2TsB(:,indx+3) ...
                         + d2zdxdyp(:,p3+[1:3])*dypdth(:,jj)*dxdri(pp);
      end

   end      

end

% d2TsB/(dr2 dth2).
% The sign on dx/dr_i is positive.  
for ii = 1:3

   i36 = 36*(ii+6-1);  % r2

   dxdri = dxdr(:,ii);
   dzdri = dzdx*dxdri;

   for jj = 1:3

      j3 = 3*(jj+9-1);  % th2
      indx = i36 + j3;

      for pp = 1:3
         p3 = 3*(pp-1);
         term1 = d2ydxdz(:,p3+[1:3])*dzdth(:,jj+3)*dxdri(pp);
         term2 = d2ydz2(:,p3+[1:3])*dzdth(:,jj+3)*dzdri(pp);
         term3 = dydz*d2zdxdyp(:,p3+[1:3])*dypdth(:,jj+3)*dxdri(pp);
         d2TsB(:,indx+2) = d2TsB(:,indx+2) + term1 + term2 + term3;
         d2TsB(:,indx+3) = d2TsB(:,indx+3) ...
                         + d2zdxdyp(:,p3+[1:3])*dypdth(:,jj+3)*dxdri(pp);
      end

   end      

end

% d2TsB/(dth1 dr1).
% The sign on dx/dr_j is negative.  
for ii = 1:3

   i36 = 36*(ii+3-1);  % th1

   for jj = 1:3

      j3 = 3*(jj-1);  % r1
      indx = i36 + j3;

      dxdrj = -dxdr(:,jj);
      dzdrj = dzdx*dxdrj;

      for pp = 1:3
         p3 = 3*(pp-1);
         term1 = d2ydxdz(:,p3+[1:3])*dzdth(:,ii)*dxdrj(pp);
         term2 = d2ydz2(:,p3+[1:3])*dzdth(:,ii)*dzdrj(pp);
         term3 = dydz*d2zdxdyp(:,p3+[1:3])*dypdth(:,ii)*dxdrj(pp);
         d2TsB(:,indx+2) = d2TsB(:,indx+2) + term1 + term2 + term3;
         d2TsB(:,indx+3) = d2TsB(:,indx+3) ...
                         + d2zdxdyp(:,p3+[1:3])*dypdth(:,ii)*dxdrj(pp);
      end

   end      

end

% d2TsB/(dth1 dr2).
% The sign on dx/dr_j is positive.  
for ii = 1:3

   i36 = 36*(ii+3-1);  % th1

   for jj = 1:3

      j3 = 3*(jj+6-1);  % r2
      indx = i36 + j3;

      dxdrj = dxdr(:,jj);
      dzdrj = dzdx*dxdrj;

      for pp = 1:3
         p3 = 3*(pp-1);
         term1 = d2ydxdz(:,p3+[1:3])*dzdth(:,ii)*dxdrj(pp);
         term2 = d2ydz2(:,p3+[1:3])*dzdth(:,ii)*dzdrj(pp);
         term3 = dydz*d2zdxdyp(:,p3+[1:3])*dypdth(:,ii)*dxdrj(pp);
         d2TsB(:,indx+2) = d2TsB(:,indx+2) + term1 + term2 + term3;
         d2TsB(:,indx+3) = d2TsB(:,indx+3) ...
                         + d2zdxdyp(:,p3+[1:3])*dypdth(:,ii)*dxdrj(pp);
      end

   end      

end

% d2TsB/(dth2 dr1).
% The sign on dx/dr_j is negative.  
for ii = 1:3

   i36 = 36*(ii+9-1);  % th2

   for jj = 1:3

      j3 = 3*(jj-1);  % r1
      indx = i36 + j3;

      dxdrj = -dxdr(:,jj);
      dzdrj = dzdx*dxdrj;

      for pp = 1:3
         p3 = 3*(pp-1);
         term1 = d2ydxdz(:,p3+[1:3])*dzdth(:,ii+3)*dxdrj(pp);
         term2 = d2ydz2(:,p3+[1:3])*dzdth(:,ii+3)*dzdrj(pp);
         term3 = dydz*d2zdxdyp(:,p3+[1:3])*dypdth(:,ii+3)*dxdrj(pp);
         d2TsB(:,indx+2) = d2TsB(:,indx+2) + term1 + term2 + term3;
         d2TsB(:,indx+3) = d2TsB(:,indx+3) ...
                         + d2zdxdyp(:,p3+[1:3])*dypdth(:,ii+3)*dxdrj(pp);
      end

   end      

end

% d2TsB/(dth2 dr2).
% The sign on dx/dr_j is positive.  
for ii = 1:3

   i36 = 36*(ii+9-1);  % th2

   for jj = 1:3

      j3 = 3*(jj+6-1);  % r2
      indx = i36 + j3;

      dxdrj = dxdr(:,jj);
      dzdrj = dzdx*dxdrj;

      for pp = 1:3
         p3 = 3*(pp-1);
         term1 = d2ydxdz(:,p3+[1:3])*dzdth(:,ii+3)*dxdrj(pp);
         term2 = d2ydz2(:,p3+[1:3])*dzdth(:,ii+3)*dzdrj(pp);
         term3 = dydz*d2zdxdyp(:,p3+[1:3])*dypdth(:,ii+3)*dxdrj(pp);
         d2TsB(:,indx+2) = d2TsB(:,indx+2) + term1 + term2 + term3;
         d2TsB(:,indx+3) = d2TsB(:,indx+3) ...
                         + d2zdxdyp(:,p3+[1:3])*dypdth(:,ii+3)*dxdrj(pp);
      end

   end      

end

% d2TsB/(dth1 dth1).
for ii = 1:3

   i36 = 36*(ii+3-1);  % th1

   for jj = 1:3

      j3 = 3*(jj+3-1);  % th1
      indx = i36 + j3;

      for pp = 1:3
         p6 = 6*(pp-1);
         p3 = 3*(pp-1);
         dip = (ii == pp);
         term1 = d2ydz2(:,p3+[1:3])*dzdth(:,jj)*dzdth(pp,ii);
         term2 = d2zdyp2(:,p3+[1:3])*dypdth(:,jj)*dypdth(pp,ii);
         term3 = dzdyp*d2ypdth2(:,p6+jj)*dip;
         d2TsB(:,indx+2) = d2TsB(:,indx+2) + term1 + dydz*term2 + dydz*term3;
         d2TsB(:,indx+3) = d2TsB(:,indx+3) + term2 + term3;
      end

   end      

end

% d2TsB/(dth1 dth2).
for ii = 1:3

   i36 = 36*(ii+3-1);  % th1

   for jj = 1:3

      j3 = 3*(jj+9-1);  % th2
      indx = i36 + j3;

      for pp = 1:3
         p6 = 6*(pp-1);
         p3 = 3*(pp-1);
         dip = (ii == pp);
         term1 = d2ydz2(:,p3+[1:3])*dzdth(:,jj+3)*dzdth(pp,ii);
         term2 = d2zdyp2(:,p3+[1:3])*dypdth(:,jj+3)*dypdth(pp,ii);
         term3 = dzdyp*d2ypdth2(:,p6+jj+3)*dip;
         d2TsB(:,indx+2) = d2TsB(:,indx+2) + term1 + dydz*term2 + dydz*term3;
         d2TsB(:,indx+3) = d2TsB(:,indx+3) + term2 + term3;
      end

   end      

end

% d2TsB/(dth2 dth1).
for ii = 1:3

   i36 = 36*(ii+9-1);  % th2

   for jj = 1:3

      j3 = 3*(jj+3-1);  % th1
      indx = i36 + j3;

      for pp = 1:3
         p6 = 6*(pp+3-1);
         p3 = 3*(pp-1);
         dip = (ii == pp);
         term1 = d2ydz2(:,p3+[1:3])*dzdth(:,jj)*dzdth(pp,ii+3);
         term2 = d2zdyp2(:,p3+[1:3])*dypdth(:,jj)*dypdth(pp,ii+3);
         term3 = dzdyp*d2ypdth2(:,p6+jj)*dip;
         d2TsB(:,indx+2) = d2TsB(:,indx+2) + term1 + dydz*term2 + dydz*term3;
         d2TsB(:,indx+3) = d2TsB(:,indx+3) + term2 + term3;
      end

   end      

end

% d2TsB/(dth2 dth2).
for ii = 1:3

   i36 = 36*(ii+9-1);  % th2

   for jj = 1:3

      j3 = 3*(jj+9-1);  % th2
      indx = i36 + j3;

      for pp = 1:3
         p6 = 6*(pp+3-1);
         p3 = 3*(pp-1);
         dip = (ii == pp);
         term1 = d2ydz2(:,p3+[1:3])*dzdth(:,jj+3)*dzdth(pp,ii+3);
         term2 = d2zdyp2(:,p3+[1:3])*dypdth(:,jj+3)*dypdth(pp,ii+3);
         term3 = dzdyp*d2ypdth2(:,p6+jj+3)*dip;
         d2TsB(:,indx+2) = d2TsB(:,indx+2) + term1 + dydz*term2 + dydz*term3;
         d2TsB(:,indx+3) = d2TsB(:,indx+3) + term2 + term3;
      end

   end      

end

