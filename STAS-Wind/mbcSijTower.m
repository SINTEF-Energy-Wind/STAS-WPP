function [SpsiR,SpsiI] = mbcSijTower (Neb,Nut,nf,df,Vinf,I,Lu,Omega, ...
                                      r,zut_g,sgr_g,Tg_r,Ty_g,psi0)
% 03.08.2017  Adapted for complex step derivatives.

smoothFlag = 1;

T0 = 1/df;
dt = T0/(4*nf);

Q = buildQijTower (Neb,Nut,nf,df,Vinf,I,Lu,Omega, ...
                   r,zut_g,sgr_g,Tg_r,Ty_g);
QijR = mbcQijTower (Neb,Nut,Q,nf,df,psi0,Omega);

clear Q;

QijI = imag(QijR);
QijR = real(QijR);

% For complex step derivatives: Qij is nominally a real array, but
% now it may be complex.  We need to separate the nominal, real
% part from the complex perturbation, and FFT the two separately.
% Then:
%     Q     FFT     real(SpsiR) imag(SpsiR) real(SpsiI) imag(SpsiI)
%   real   real          X
%   real   imag                                  X
%   imag   real                      X
%   imag   imag                                              X
SpsiR = zeros(size(QijR),'single');
SpsiI = zeros(size(QijR),'single');
M = Neb;
N = size(QijR,1);
dcol = size(QijR,2)/M;
for j = 1:M
   icol = dcol*(j-1);

   % Extract the part of Qij we need.
   Qij = QijR(:,icol+1:icol+dcol);

   Store = dt*fft(Qij);

   if (smoothFlag == 1)
      Ssmooth = zeros(N,dcol);
      Ssmooth(1,1:dcol) = 0.25*Store(N,:) ...
                        + 0.50*Store(1,:) ...
                        + 0.25*Store(2,:);
      Ssmooth(2:N-1,1:dcol) = 0.25*Store(1:N-2,:) ...
                            + 0.50*Store(2:N-1,:) ...
                            + 0.25*Store(3:N  ,:);
      Ssmooth(N,1:dcol) = 0.25*Store(N-1,:) ...
                        + 0.50*Store( N ,:) ...
                        + 0.25*Store( 1 ,:);

      SpsiR(:,icol+1:icol+dcol) = real(Ssmooth);
      SpsiI(:,icol+1:icol+dcol) = imag(Ssmooth);
   end

end

clear QijR;

for j = 1:M
   icol = dcol*(j-1);

   Qij = QijI(:,icol+1:icol+dcol);

   Store = dt*fft(Qij);

   if (smoothFlag == 1)
      Ssmooth = zeros(N,dcol);
      Ssmooth(1,1:dcol) = 0.25*Store(N,:) ...
                        + 0.50*Store(1,:) ...
                        + 0.25*Store(2,:);
      Ssmooth(2:N-1,1:dcol) = 0.25*Store(1:N-2,:) ...
                            + 0.50*Store(2:N-1,:) ...
                            + 0.25*Store(3:N  ,:);
      Ssmooth(N,1:dcol) = 0.25*Store(N-1,:) ...
                        + 0.50*Store( N ,:) ...
                        + 0.25*Store( 1 ,:);

      SpsiR(:,icol+1:icol+dcol) = SpsiR(:,icol+1:icol+dcol) ...
                                + i*real(Ssmooth);
      SpsiI(:,icol+1:icol+dcol) = SpsiI(:,icol+1:icol+dcol) ...
                                + i*imag(Ssmooth);
   end

end

