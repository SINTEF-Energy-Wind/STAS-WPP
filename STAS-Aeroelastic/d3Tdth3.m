function d3T = d3Tdth3 (th,T,dT,d2T)
%
% Compute through the third derivatives of a transform matrix in terms
% of the exponential map rotations.
%
% Version:        Changes:
% --------        -------------
% 03.11.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 03.11.2017      Verified by finite difference.
%

d3T   = zeros(3,81);

normth = sqrt((th.')*th);

if (abs(normth) < eps^0.25)  % Complex magnitude is intentionally included
                             % in this norm.
   
   TTH   = vecToSpin (th);
   normth2 = normth^2;

   c = (-(1/3) + normth2/30);

   for jj = 1:3

      jc27 = 27*(jj-1);

      ej = zeros(3,1);
      ej(jj) = 1;
      Sigj = vecToSpin (ej);
      thj = th(jj);

      for kk = 1:3

         kc9 = 9*(kk-1);

         ek = zeros(3,1);
         ek(kk) = 1;
         Sigk = vecToSpin (ek);
         djk = (jj == kk);
         thk = th(kk);

         for pp = 1:3

            pc3 = 3*(pp-1);

            ep = zeros(3,1);
            ep(pp) = 1;
            Sigp = vecToSpin (ep);
            djp = (jj == pp);
            dkp = (kk == pp);
            thp = th(pp);
            
            d3T(:,jc27+kc9+pc3+[1:3]) =                          ...
                 (thk*djp + thj*dkp + thp*djk)*TTH/15            ...
               + (thj*thk*Sigp + thp*thj*Sigk + thp*thk*Sigj)/15 ...
               + c*(djk*Sigp + djp*Sigk + dkp*Sigj)              ...
               - (djk*(Sigp*TTH + TTH*Sigp)                      ...
               +  djp*(Sigk*TTH + TTH*Sigk)                      ...
               +  dkp*(Sigj*TTH + TTH*Sigj))/12                  ...
               - (thj*(Sigk*Sigp + Sigp*Sigk)                    ...
               +  thk*(Sigj*Sigp + Sigp*Sigj)                    ...
               +  thp*(Sigk*Sigj + Sigj*Sigk))/12;

         end
      end
   end

else

   TTH   = vecToSpin (th);
   dTTH  = [vecToSpin([1;0;0]) vecToSpin([0;1;0]) vecToSpin([0;0;1])];
   IT    = eye(3) - T;
   TTHIT = TTH*IT;

   magth = (th.')*th;

   for ii = 1:3

      ii3 = 3*(ii-1);
      ii27 = 27*(ii-1);

      eei = zeros(3,1);
      eei(ii) = 1;

      R = (th(ii)*TTH + vecToSpin(TTHIT(:,ii)))/magth;

      for jj = 1:3

         jj3 = 3*(jj-1);
         jj9 = 9*(jj-1);

         if (ii == jj)
            dij = 1;
         else
            dij = 0;
         end

         dRj = ((dij*TTH + th(ii)*dTTH(:,jj3+[1:3])                         ...
             +   vecToSpin((dTTH(:,jj3+[1:3])*IT-TTH*dT(:,jj3+[1:3]))*eei)) ...
             - (2*th(jj)*R))/magth;

         for kk = 1:3

            kk3 = 3*(kk-1);

            if (jj == kk)
               djk = 1;
            else
               djk = 0;
            end

            if (ii == kk)
               dik = 1;
            else
               dik = 0;
            end

            dRk = ((dik*TTH + th(ii)*dTTH(:,kk3+[1:3])                         ...
                +   vecToSpin((dTTH(:,kk3+[1:3])*IT-TTH*dT(:,kk3+[1:3]))*eei)) ...
                - (2*th(kk)*R))/magth;

            d2R = ((-2*djk*R - 2*th(kk)*dRj - 2*th(jj)*dRk)          ...
                +  (dij*dTTH(:,kk3+[1:3]) + dik*dTTH(:,jj3+[1:3])    ...
                +   vecToSpin(-dTTH(:,jj3+[1:3])*dT(:,kk3+[1:3])*eei ...
                -             dTTH(:,kk3+[1:3])*dT(:,jj3+[1:3])*eei  ...
                -             TTH*d2T(:,jj9+kk3+[1:3])*eei)))        ...
                / magth;

            d3T(:,ii27+jj9+kk3+[1:3]) = d2R*T + dRj*dT(:,kk3+[1:3]) ...
                                      + dRk*dT(:,jj3+[1:3]) + R*d2T(:,jj9+kk3+[1:3]);

         end

      end

   end

end