function [A,B,C,Phi,slap] = Ytransform (aa,bb,cc)

Nx = size(aa,1);

[slap,shp,ifrq] = eigVal_silent (aa);

Phi = shp(:,ifrq);
Psi = inv(Phi);

ii = [1:Nx].';
jj = [1:Nx].';
ss = slap;
Lam = sparse(ii,jj,ss,Nx,Nx);

PB = Psi*bb;
CP = cc*Phi;

iYY = speye(Nx);
Nosc = 0;
for im1 = 1:floor(Nx/2)

   imn = Nx - (im1-1);

   if (abs(imag(slap(im1))) > 0)
      Nosc = Nosc + 1;
      Phi(:,imn) = conj(Phi(:,im1));
      Psi(imn,:) = conj(Psi(im1,:));
      slap(imn) = conj(slap(im1));
      iYY([im1 imn],[im1 imn]) = [1 1;-i i];
   end

end
Nexp = Nx - 2*Nosc;
Nmds = Nosc + Nexp;

PB = Psi*bb;
CP = cc*Phi;

YY = inv(iYY);
A = real(iYY*Lam*YY);
B = real(iYY*PB);
C = real(CP*YY);