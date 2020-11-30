
clear;

load 'stress.txt';
load 'csb.txt';

WW = 0.868860523;

freqs = stress(:,1);
Nf = size(freqs,1);
df = freqs(2) - freqs(1);
dt = 1/(Nf*df);
T = (Nf-1)*dt;

Nrep = 200;  % Number of repititions.

Sf = stress(:,6) + i*stress(:,7);
Se = stress(:,8) + i*stress(:,9);
cf = csb(:,2) + i*csb(:,3);
ce = csb(:,4) + i*csb(:,5);

Sfp = cf.*conj(cf)/df;
Sep = ce.*conj(ce)/df;

% Isolate parts of the spectra.
P = WW/(2*pi);
w = round(P/df);
w2 = round(w/2);

indL = Nf/2-w2+[1:2*w2+1].';
SfL = Sf(indL);
SeL = Se(indL);

cfL = cf(indL);
ceL = ce(indL);

indS = [Nf/2-w2-w+[1:w] Nf/2+1+w2+[1:w]].';
SfS = Sf(indS);
SeS = Se(indS);

cfP = cf(indS);
ceP = ce(indS);

indH = [[1:Nf/2-w2-w] [Nf/2+2+w2+w:Nf]].';
SfH = (Sf(indH) + Sfp(indH));
SeH = (Se(indH) + Sep(indH));

% Output some time series.
fidfL   = fopen('data_fL.txt','w');
fideL   = fopen('data_eL.txt','w');
fidfS   = fopen('data_fS.txt','w');
fideS   = fopen('data_eS.txt','w');
fidfP   = fopen('data_fP.txt','w');
fideP   = fopen('data_eP.txt','w');
fidfH   = fopen('data_fH.txt','w');
fideH   = fopen('data_eH.txt','w');

fidfSP  = fopen('data_fSP.txt','w');
fideSP  = fopen('data_eSP.txt','w');
fidfSPH = fopen('data_fSPH.txt','w');
fideSPH = fopen('data_eSPH.txt','w');
fidf    = fopen('data_f.txt','w');
fide    = fopen('data_e.txt','w');

for irep = 1:Nrep

printf('%5d\n',irep);
fflush(stdout);

   a1 = -pi + 2*pi*rand(Nf/2+1,1);
   a2 = -pi + 2*pi*rand(Nf/2+1,1);
   a3 = -pi + 2*pi*rand(Nf/2+1,1);
   a4 = -pi + 2*pi*rand(Nf/2+1,1);
   a5 = -pi + 2*pi*rand(Nf/2+1,1);
   a6 = -pi + 2*pi*rand(Nf/2+1,1);
   alf1 = [a1; -a1(Nf/2:-1:2)];
   alf2 = [a2; -a2(Nf/2:-1:2)];
   alf3 = [a3; -a3(Nf/2:-1:2)];
   alf4 = [a4; -a4(Nf/2:-1:2)];
   alf5 = [a5; -a5(Nf/2:-1:2)];
   alf6 = [a6; -a6(Nf/2:-1:2)];
   for it = 1:Nf

      t = (it-1)*dt;
   
      ftL = real(sum(sqrt(abs(SfL)*df).*exp(i*(2*pi*freqs(indL)*t + alf1(indL)))));
      fprintf(fidfL,'%+5.6e %+5.6e\n',t,ftL);

      etL = real(sum(sqrt(abs(SeL)*df).*exp(i*(2*pi*freqs(indL)*t + alf2(indL)))));
      fprintf(fideL,'%+5.6e %+5.6e\n',t,etL);

      ftS = real(sum(sqrt(abs(SfS)*df).*exp(i*(2*pi*freqs(indS)*t + alf3(indS)))));
%      ftS2 = sum(sqrt(abs(SfS)*df).*sin(2*pi*freqs(indS)*t + alf3(indS)));
      fprintf(fidfS,'%+5.6e %+5.6e\n',t,ftS);

      etS = real(sum(sqrt(abs(SeS)*df).*exp(i*(2*pi*freqs(indS)*t + alf4(indS)))));
      fprintf(fideS,'%+5.6e %+5.6e\n',t,etS);

      ftH = real(sum(sqrt(abs(SfH)*df).*exp(i*(2*pi*freqs(indH)*t + alf5(indH)))));
      fprintf(fidfH,'%+5.6e %+5.6e\n',t,ftH);

      etH = real(sum(sqrt(abs(SeH)*df).*exp(i*(2*pi*freqs(indH)*t + alf6(indH)))));
      fprintf(fideH,'%+5.6e %+5.6e\n',t,etH);   

      ftP = sum(cfP.*exp(i*2*pi*freqs(indS)*t));
      fprintf(fidfP,'%+5.6e %+5.6e\n',t,ftP);   

      etP = sum(ceP.*exp(i*2*pi*freqs(indS)*t));
      fprintf(fideP,'%+5.6e %+5.6e\n',t,etP);

      ftSP = ftS + ftP;
      fprintf(fidfSP,'%+5.6e %+5.6e\n',t,ftSP);

      etSP = etS + etP;
      fprintf(fideSP,'%+5.6e %+5.6e\n',t,etSP);

      ftSPH = ftS + ftP + ftH;
      fprintf(fidfSPH,'%+5.6e %+5.6e\n',t,ftSPH);

      etSPH = etS + etP + etH;
      fprintf(fideSPH,'%+5.6e %+5.6e\n',t,etSPH);

      ft = ftL + ftS + ftP + ftH;
      fprintf(fidf,'%+5.6e %+5.6e\n',t,ft);

      et = etL + etS + etP + etH;
      fprintf(fide,'%+5.6e %+5.6e\n',t,et);

   end % it

end % irep
fclose('all');

nm1  = '_f';
nm2  = '_fL';
nm3  = '_fS';
nm4  = '_fP';
nm5  = '_fH';
nm6  = '_fSP';
nm7  = '_fSPH';
nm8  = '_e';
nm9  = '_eL';
nm10 = '_eS';
nm11 = '_eP';
nm12 = '_eH';
nm13 = '_eSP';
nm14 = '_eSPH';

rs1  = '_30_0100';
rs2  = '_20_0050';
rs3  = '_20_0050';
rs4  = '_05_0020';
rs5  = '_20_0050';
rs6  = '_20_0050';
rs7  = '_20_0050';
rs8  = '_30_0100';
rs9  = '_02_0005';
rs10 = '_02_0005';
rs11 = '_20_0050';
rs12 = '_20_0050';
rs13 = '_20_0050';
rs14 = '_30_0100';

for jnm = 1:14

   eval(['nm = nm' int2str(jnm) ';']);
printf('%s\n',nm);
fflush(stdout);
   eval(['load data' nm '.txt']);
   eval(['dat = data' nm ';']);

   eval(['rs = rs' int2str(jnm) ';']);

   for irep = 1:Nrep

      irow = Nf*(irep-1);

      dd = dat(irow+[1:Nf],:);

      fid = fopen('data.txt','w');
      for it = 1:Nf
         fprintf(fid,'%+5.6e %+5.6e\n',dd(it,1),dd(it,2));
      end
      fclose(fid);

      eval(["system('rainflow" rs ".exe');"]);

      load 'ncum.txt';
      load 'ncyc.txt';

      if (irep == 1)
         nbin = size(ncyc,1);
         sbin = ncyc(nbin:-1:1,1);
         ds = sbin(2) - sbin(1);
         ntot = zeros(nbin,3);
         ntot(:,1) = sbin;
      end

      cyc = ncyc(nbin:-1:1,2);
      cum = ncum(nbin:-1:1,2);
      ntot(:,2) = ntot(:,2) + cyc/Nrep;
      ntot(:,3) = ntot(:,3) + cum/Nrep;

   end

   eval(["fid = fopen('ncyc" nm ".txt','w');"]);
   for ib = 1:nbin
      fprintf(fid,'%+5.6e %+5.6e %+5.6e\n', ...
              0.5*ntot(ib,1),ntot(ib,2)/(T*(0.5*ds)),ntot(ib,3)/(T*(0.5*ds)));
   end
   fclose(fid);

end




