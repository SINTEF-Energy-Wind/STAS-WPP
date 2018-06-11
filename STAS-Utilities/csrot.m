function T = csrot (axan)
%
% axan is axis, angle pairs: [axis1 angle1;axis2 angle2;...]
%
% axis: x=1, y=2, z=3.  Use positive values for premultiplication,
% (rotate about original axis) negative values for postmultiplication
% (rotate about updated axis).
%
% The function returns T_2^1.  T_1^2 = (T_2^1)^T.
%

Nop = size(axan,1);

T = eye(3);
n = 0;
while (n < Nop) 

   n++;  

   ax =  abs(axan(n,1));
   cs = sign(axan(n,1));

   ca =  cos(axan(n,2));
   sa =  sin(axan(n,2));

   if (ax == 1)
      TT = [1   0   0;  ...
            0  ca -sa;  ...
            0  sa  ca];
   elseif (ax == 2)
      TT = [ca  0  sa;  ...
             0  1   0;  ...
           -sa  0  ca];
   elseif (ax == 3)
      TT = [ca -sa  0;  ...
            sa  ca  0;  ...
             0   0  1];
   end

   if (cs >= 0)
      T = TT*T;
   else
      T = T*TT;
   end

endwhile

