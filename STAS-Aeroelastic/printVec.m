function Vc = printVec(V)

N = length(V);
Nr = N/6;

Vc = zeros(Nr,7);

for i = 1:Nr
   idof = 6*(i-1);
   Vc(i,:) = [i V(idof+1) V(idof+2) V(idof+3) V(idof+4) V(idof+5) V(idof+6)];
end