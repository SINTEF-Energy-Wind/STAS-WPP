function PnFlex (xng)

Nnod = size(xng,1)/3;

fid = fopen('Pnod.txt','w');

for inod = 1:Nnod

   ic3 = 3*(inod-1);

   fprintf(fid,'%+5.6e %+5.6e %+5.6e\n',xng(ic3+1),xng(ic3+2),xng(ic3+3));

end

fclose('all');