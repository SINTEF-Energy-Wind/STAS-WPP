function printMatrix (fid,M)

Nr = size(M,1);
Nc = size(M,2);

for ir = 1:Nr
   fprintf(fid,'%+5.6e',M(ir,1));
   for ic = 2:Nc
      fprintf(fid,' %+5.6e',M(ir,ic));
   end
   fprintf(fid,'\n');
end
