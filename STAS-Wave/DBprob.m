fid = fopen('prob.txt','w');
fid2 = fopen('table.txt','w');

dTp = 0.1;
dHs = 0.1;
%Tp = [1:dTp:24]';
%Hs = [0.1:dHs:20]';

dTp = 1;
dHs = 0.5;
Tp = [2:dTp:20]';
Hs = [0.5:dHs:10]';


Nh = size(Hs,1);
Nt = size(Tp,1);

ph = zeros(Nh,1);
p = zeros(Nt,Nh);

beta = 1.315;
rho = 1.546;
eta = 2.802;
alpha = 0.591;
theta = 0.313;
a1 = 0.889;
a2 = 0.913;
a3 = 0.300;
b1 = 0.005;
b2 = 0.123;
b3 = 0.486;
mu = a1 + a2*(Hs.^a3);
sigma = sqrt(b1 + b2*exp(-b3*Hs));

hh = Hs(Hs <= eta);
ph(Hs<=eta) = exp(-((log(hh) - theta).^2)/(2*(alpha^2))) ...
          ./ (sqrt(2*pi)*alpha*hh);

hh = Hs(Hs > eta);
ph(Hs>eta) = (beta/rho)*((hh/rho).^(beta-1)).*exp(-(hh/rho).^beta);

for ih = 1:Nh

   p(:,ih) = ph(ih) ...
           * exp(-((log(Tp) - mu(ih)).^2)./(2*(sigma(ih).^2))) ...
          ./ (sqrt(2*pi)*sigma(ih).*Tp);

end

p(abs(p)<1e-12) = 0;

sump = sum(sum(p)*dTp*dHs);

for ih = 1:Nh
   if (ih == 1)
      fprintf(fid2,'%5.6e',0);
      for it = 1:Nt
         fprintf(fid2,' %+5.6e',Tp(it));
      end
      fprintf(fid2,'\n');
   end
   fprintf(fid2,'%+5.6e',Hs(ih));
   for it = 1:Nt
%      fprintf(fid,'%+5.6e %+5.6e %+5.6e\n', ...
%              Tp(it),Hs(ih),p(it,ih));
%      fprintf(fid,'%+5.6e %+5.6e %+5.6e\n', ...
%              Tp(it),Hs(ih),max(log10(p(it,ih)),-12));
      fprintf(fid,'%+5.6e %+5.6e %+5.6e\n', ...
              Tp(it),Hs(ih),p(it,ih)*dTp*dHs/sump);
      fprintf(fid2,' %+5.6e',p(it,ih)*dTp*dHs/sump);

   end
   fprintf(fid2,'\n');
end

sump

fclose('all');