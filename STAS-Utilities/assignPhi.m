function Phi = assignPhi (jc,wc,Np);

Nj1 = size(jc,1);
Nnz = Nj1*size(jc,2);
ir = reshape (jc,Nnz,1);
ind = [0:Nj1:Nnz-Nj1+1].';
ic = zeros(Nnz,1);
for irep = 1:Nj1
   ic(irep+ind) = [1:size(jc,2)].';
end
ss = reshape (wc,Nnz,1);
Phi = sparse(ir,ic,ss,Np,size(jc,2));