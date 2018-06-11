function pdot = RiccatiQS (p,t)



% CAUTION, NOT ADAPTED FOR COMPLEX STEP DERIVATIVES.



global AA_Ric BC_Ric QQ_Ric RR_Ric

N = sqrt(size(p,1));

Pmat = zeros(N,N);
for jj = 1:N
   ind = N*(jj-1);
   Pmat(:,jj) = p(ind+[1:N]);
end

Pmat = 0.5*(Pmat + Pmat.');
RinvBCP = balanceDiv(RR_Ric,(BC_Ric.')*Pmat);

Pmdot = Pmat*AA_Ric + (AA_Ric.')*Pmat - Pmat*BC_Ric*RinvBCP + QQ_Ric;

% Real: prevent small imaginary residuals from triggering warnings.
pdot = zeros(N^2,1);
for jj = 1:N
   ind = N*(jj-1);
   pdot(ind+[1:N]) = real(Pmdot(:,jj));
end
