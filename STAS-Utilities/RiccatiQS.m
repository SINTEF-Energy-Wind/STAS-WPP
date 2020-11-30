function [Pdot,dRic] = RiccatiQS (linFlag,Pin,A,B,Q,R)

N = size(Pin,1);

P = 0.5*(Pin + Pin.');
AT = A.';
BRiB = B*balanceDiv(R,(B.'));
Pdot = P*A + AT*P - P*BRiB*P + Q;
BRiBP = BRiB*P;
PBRiB = P*BRiB;

dRic = sparse(N^2,N^2);
if (linFlag == 1)

   for ii = 1:N
      irn = N*(ii-1);
      for jj = 1:N
         jcn = N*(jj-1);
         icol = jj + N*(ii-1);
         mat = sparse(N,N);
         mat(jj,:) = A(jj,:) - BRiBP(jj,:);
         mat(:,ii) = mat(:,ii) + AT(:,ii) - PBRiB(:,ii);
         dRic(:,icol) = reshape(mat,N^2,1);
      end
   end   

end

