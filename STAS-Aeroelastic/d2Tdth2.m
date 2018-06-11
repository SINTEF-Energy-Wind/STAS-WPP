function d2T = d2Tdth2 (th,T,dT)
%
% Compute the second derivatives of a transform matrix in terms of the
% exponential map rotations.
%
%                            --                            --                            -- --
%    d2T             thj     |                             |                       dT     |  |
% --------- = -2 ----------- | dij*TH + thi*spin(ej) + spin| spin(ej)(I-T)ei - TH ---- ei |  |
% dthi*dthj      (th^T*th)^2 |                             |                      dthj    |  |
%                            --                            --                            -- --
%
% A special formula based on a series expansion for sin,cos is used
% for small angle inputs, which works with complex step.  The baseline
% formula suffers from underflow with complex step.
%
% Version:        Changes:
% --------        -------------
% 16.10.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 16.10.2017      Verified by complex step and finite difference.
%

d2T = zeros(3,9*3);

normth = sqrt((th.')*th);

TTH = vecToSpin(th);

if (abs(normth) < eps^0.25)  % Complex magnitude is intentionally included
                             % in this norm.

   normth2 = normth^2;

   for ii = 1:3

      ei = zeros(3,1);
      ei(ii) = 1;
      spi = vecToSpin(ei);
      sTTsi = spi*TTH + TTH*spi;

      for jj = 1:3

         ic3 = 3*(jj-1);
         ic = 9*(ii-1) + ic3;

         ej = zeros(3,1);
         ej(jj) = 1;
         spj = vecToSpin(ej);

         dij = (ii==jj);  % 1 if so, 0 if not.

         d2T(:,ic+[1:3]) = -dij*TTH/3 - (th(jj)*spi + th(ii)*spj)/3        ...
                         - dij*TTH*TTH/12                                  ...
                         - th(jj)*sTTsi/12 - th(ii)*(spj*TTH + TTH*spj)/12 ...
                         + (0.5 - normth2/24)*(spj*spi + spi*spj); 

      end

   end

else

   for ii = 1:3

      ei = zeros(3,1);
      ei(ii) = 1;

      term1 = th(ii)*TTH;
      term2 = vecToSpin(TTH*(eye(3) - T)*ei);
      mat = ((term1 + term2)/(normth^2));

      for jj = 1:3

         ic3 = 3*(jj-1);
         ic = 9*(ii-1) + ic3;

         ej = zeros(3,1);
         ej(jj) = 1;
         spj = vecToSpin(ej);

         dd = (ii==jj);  % 1 if so, 0 if not.
   
         vec = spj*(eye(3) - T)*ei - TTH*dT(:,ic3+[1:3])*ei;

         d2T(:,ic+[1:3]) = -2*(th(jj)/((th.'*th)^2))*(term1 + term2)*T ...
                         + (dd*TTH + th(ii)*spj + vecToSpin(vec))*T/(normth^2) ...
                         + mat*dT(:,ic3+[1:3]);

      end

   end

end