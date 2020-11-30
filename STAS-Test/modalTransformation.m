function [Am,Bm,Cm] = modalTransformation (A,B,C,freqs,tol,Yflg,wrtflg)
%
% Compute the eigenmodes of the matrix A, and retain the number
% that contribute to the transfer functions with relative magnitude
% "tol" or more, at the specified frequencies "freqs".
%
% Version:        Changes:
% --------        -------------
% 07.07.2020      Original code.
%
% Version:        Verification:
% --------        -------------
% 07.07.2020      
%
% Inputs:
% -------
% A, B, C         : State matrices.
% freqs           : Vector of frequencies at which to evaluate 
%                   TFs.
% tol             : Tolerance for TF errors.
% Yflg            : = 0: Y=I, output matrices may be complex.
%                   = 1: Perform a Y transform for real outputs.
% wrtflg          : = 0: Do not write anything.
%                   = 1: Write ranked mode table. 
%
% Outputs:
% --------
% Am, Bm, Cm      : Modal matrices.

Nx = size(A,1);
Nu = size(B,2);
Ny = size(C,1);
Nf = size(freqs,1);
NTFs = Nu*Ny;

[slap,shp,ifrq] = eigVal_silent (A);

% Diagonalize.
Phim = shp(:,ifrq);
Psim = inv(Phim);

% Some of the very closely spaced eigenfrequencies may lie out of order.
% This causes problems when I assume that conjugate modes are mirrored
% in the sorted mode shape matrix.  The solution is to force complex
% conjugacy by copying half the complex mode shapes and eigenvalues to
% the other half.  This is done in the 'for' loop below.
iY = speye(Nx);
Nosc = 0;
for im1 = 1:floor(Nx/2)

   imn = Nx - (im1-1);

   if (abs(imag(slap(im1))) > 0)
      Nosc = Nosc + 1;
      Phim(:,imn) = conj(Phim(:,im1));
      Psim(imn,:) = conj(Psim(im1,:));
      slap(imn) = conj(slap(im1));
      if (Yflg == 1)
         iY([im1 imn],[im1 imn]) = [1 1;-i i];
      end
   end

end
Nexp = Nx - 2*Nosc;
Nmds = Nosc + Nexp;

ii = [1:Nx].';
jj = ii;
PAP = sparse (ii,jj,slap,Nx,Nx);  % = Psim*A*Phim.
PB = Psim*B;
CP = C*Phim;

% Evaluate the mode-by-mode TFs at the specified frequencies.
for ifreq = 1:Nf

   f = freqs(ifreq);
   w = 2*pi*f;

   % Disturbances-to-states.
   iwslap = i*w - slap.';

   for jinp = 1:Nu
      for jout = 1:Ny
         icn = Nx*(jout-1) + Nx*Ny*(jinp-1);
         Hmod(ifreq,icn+[1:Nx]) = CP(jout,:).*(PB(:,jinp).')./iwslap;
      end
   end

end

H = zeros (Nf,NTFs);
for itf = 1:NTFs
   icn = Nx*(itf-1);
   H(:,itf) = sum(Hmod(:,icn+[1:Nx]),2);
end

Hmag = abs(H);
Hang = atan2(imag(H),real(H))/pi;

% Compute projections.  Use the max error in the projection over freqs
% as a metric.
proj = zeros(Nf,NTFs*Nmds);
for itf = 1:NTFs

   icr = Nx*(itf-1);
   icn = Nmds*(itf-1);
   hvec = H(:,itf)./(conj(H(:,itf)).*H(:,itf));

   proj(:,icn+[1:Nosc]) = (conj(Hmod(:,icr+[1:Nosc]))             ...
                        +  conj(Hmod(:,icr+[Nx:-1:Nx-Nosc+1]))).*hvec;
   proj(:,icn+Nosc+[1:Nexp]) = conj(Hmod(:,icr+Nosc+[1:Nexp])).*hvec;

end

pval = zeros(Nf,NTFs*Nmds);
sindx = zeros(Nf,NTFs*Nmds);
for itf = 1:NTFs
   icn = Nmds*(itf-1);
% Projections aligned with the full TF.
%   [pval(:,icn+[1:Nmds]),sindx(:,icn+[1:Nmds])] = ...
%                       sort (abs(real(proj(:,icn+[1:Nmds]))),2,'descend');
% Magnitude, both aligned and orthogonal components.
   [pval(:,icn+[1:Nmds]),sindx(:,icn+[1:Nmds])] = ...
                       sort (abs(proj(:,icn+[1:Nmds])),2,'descend');

   for ifreq = 1:Nf
      pval(ifreq,icn+[1:Nmds]) = proj(ifreq,icn+sindx(ifreq,icn+[1:Nmds]));
   end
end

pcum = zeros(Nf,NTFs*Nmds);
for itf = 1:NTFs

   icn = Nmds*(itf-1);

   pcum(:,icn+1) = pval(:,icn+1);
   for imd = 2:Nmds
      pcum(:,icn+imd) = pcum(:,icn+imd-1) + pval(:,icn+imd);
   end

end

mlist = sindx(abs(pval) >= tol);
[mls, inds] = sort(mlist);
mls = mls(mls > 0);
Nmls = size(mls,1);

ml = zeros(Nmls,1);
kk = 0;
prev = 0;
for jj = 1:Nmls
   if (mls(jj) ~= prev)
      kk = kk + 1;
      ml(kk) = mls(jj);
      prev = ml(kk);
   end
end

ml = ml(ml > 0);
Nml = size(ml,1);

% Get the full set of indices (+ and - frequencies).
mlf = zeros(2*Nml,1);
kk = 0;
for jj = 1:Nml
   if (ml(jj) <= Nosc)
      kk = kk + 1;
      mlf(kk) = ml(jj);
      kk = kk + 1;
      mlf(kk) = Nx - ml(jj) + 1;
   else
      kk = kk + 1;
      mlf(kk) = ml(jj);
   end
end

[mlf,indm] = sort (mlf);
mlf = mlf(mlf > 0);
Nxi = size(mlf,1);

Am = sparse (diag(slap(mlf)));
Bm = PB(mlf,:);
Cm = CP(:,mlf);

if (Yflg == 1)

   % Perform the Y transform to get real matrices.
   iYY = iY(mlf,mlf);
   YY = inv(iYY);
   Am = real(iYY*Am*YY);
   Bm = real(iYY*Bm);
   Cm = real(Cm*YY);

end

