clear;

dtW = 5;
NtW = 2000;

Nturb = 32;

W = zeros(NtW,2*Nturb);
for iturb = 1:Nturb

   ic2 = 2*(iturb-1);

   txt = ["./Wind/10mps_NTM_TotalControlWPP_VS_T" int2str(iturb) ".hh"];
   dat = importdata (txt,' ',8);
   Vmag = dat.data(:,2);
   thdeg = dat.data(:,3);
   W(:,ic2+1) = Vmag;
   W(:,ic2+2) = thdeg*pi/180;

end

save('-binary','WindTS.bin','W');