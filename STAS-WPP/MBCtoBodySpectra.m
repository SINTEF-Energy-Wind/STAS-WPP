function SB = MBCtoBodySpectra (opt,Spsi,azi,W,df)
%
% Convert an input spectrum in multi-blade or nacelle coordinates to
% an output spectrum in rotating coordinates.  Options are provided
% for transforming from multi-blade coordinates to each blade, and
% from nacelle to driveshaft coordinates.
%
% The input and output spectra are packed such that each column is
% a spectrum, the rows being the frequencies, in the order
% 0 ... N/2-1, -N/2, ... -1.
% The columns are arranged in the order
% S11, S12, S13, S21, S22, S23, S31, S32, S33.
% If an MBC transform is being performed then the indices represent
% the blade number, and if a driveshaft transform is being performed
% then the indices represent 1=x, 2=y, 3=z.
%
% Version:        Changes:
% --------        -------------
% 22.08.2020      Original code.
%
% Version:        Verification:
% --------        -------------
% 22.08.2020      Verified against an independent calculation using
%                 wind fields from buildSij.
%
% Inputs:
% -------
% opt             : = 1 for MBC to body,
%                   = 2 for driveshaft.
% Spsi            : the MBC spectral matrix, packed as indicated above.
% azi             : azimuth angle at t = 0 (usually 0).
% W               : Rotor speed, rad/s.
% df              : Frequency bin, Hz.
%
% Outputs:
% --------
% SB              : the transformed spectra, packed like Spsi.

Nf = size(Spsi,1);
ww = W/(2*pi*df);

ei0 = exp(i*azi);

if (opt == 1)  % MBC.
   ei2  = exp(i*2*pi/3);
   ei4  = exp(i*4*pi/3);
   E0 = [1 0 0; ...
         1 0 0; ...
         1 0 0];
   E = 0.5*[0,   1,     -i; ...
            0, ei2, -i*ei2; ...
            0, ei4, -i*ei4]*ei0;
elseif (opt == 2)  % Driveshaft.
   E0 = [0 0 0; ...
         0 0 0; ...
         0 0 1];
   E = 0.5*[1, -i,  0; ...
            i,  1,  0; ...
            0,  0,  0]*ei0;
elseif (opt == -1)  % Inverse MBC.
   ei2  = exp(i*2*pi/3);
   ei4  = exp(i*4*pi/3);
   third = 1/3;
   E0 = third*[1 1 1; ...
               0 0 0; ...
               0 0 0];
   E = third*[0,      0,      0; ...
              1,    ei2,    ei4; ...
             -i, -i*ei2, -i*ei4];
elseif (opt == -2)  % Inverse driveshaft.
   E0 = [0 0 0; ...
         0 0 0; ...
         0 0 1];
   E = 0.5*[1,  i,  0; ...
           -i,  1,  0; ...
            0,  0,  0]*ei0;
end

Estar = conj(E);
T0 = E0 + E + Estar;

for ir = 1:3
   for ic = 1:3
      icol = 3*(ir-1) + ic;
      SB(:,icol) = T0(ir,1)*Spsi(:,1)*E0(ic,1) ...
                 + T0(ir,1)*Spsi(:,2)*E0(ic,2) ...
                 + T0(ir,1)*Spsi(:,3)*E0(ic,3) ...
                 + T0(ir,2)*Spsi(:,4)*E0(ic,1) ...
                 + T0(ir,2)*Spsi(:,5)*E0(ic,2) ...
                 + T0(ir,2)*Spsi(:,6)*E0(ic,3) ...
                 + T0(ir,3)*Spsi(:,7)*E0(ic,1) ...
                 + T0(ir,3)*Spsi(:,8)*E0(ic,2) ...
                 + T0(ir,3)*Spsi(:,9)*E0(ic,3);
   end
end

% Re-map from one frequency to another, including interpolation
% associated with non-integer w.
val = ww - floor(ww);
if (val == 0.)
   alf = 0.;
else
   alf = 1 - val;
end
ifreqs = [1:Nf].';
ilow = floor(ifreqs - ww);
ihi = ilow + 1;
ind = (ilow <= 0);
ilow(ind) = ilow(ind) + Nf;
ind = (ihi <= 0);
ihi(ind) = ihi(ind) + Nf;
ind = (ihi >= Nf+1);
ihi(ind) = ihi(ind) - Nf;
for ir = 1:3
   for ic = 1:3
      icol = 3*(ir-1) + ic;
      SB(:,icol) = SB(:,icol)                                                ...
                 + T0(ir,1)*((1-alf)*Spsi(ilow,1) + alf*Spsi(ihi,1))*E(ic,1) ...
                 + T0(ir,1)*((1-alf)*Spsi(ilow,2) + alf*Spsi(ihi,2))*E(ic,2) ...
                 + T0(ir,1)*((1-alf)*Spsi(ilow,3) + alf*Spsi(ihi,3))*E(ic,3) ...
                 + T0(ir,2)*((1-alf)*Spsi(ilow,4) + alf*Spsi(ihi,4))*E(ic,1) ...
                 + T0(ir,2)*((1-alf)*Spsi(ilow,5) + alf*Spsi(ihi,5))*E(ic,2) ...
                 + T0(ir,2)*((1-alf)*Spsi(ilow,6) + alf*Spsi(ihi,6))*E(ic,3) ...
                 + T0(ir,3)*((1-alf)*Spsi(ilow,7) + alf*Spsi(ihi,7))*E(ic,1) ...
                 + T0(ir,3)*((1-alf)*Spsi(ilow,8) + alf*Spsi(ihi,8))*E(ic,2) ...
                 + T0(ir,3)*((1-alf)*Spsi(ilow,9) + alf*Spsi(ihi,9))*E(ic,3);
   end
end

alf = ww - floor(ww);
ilow = floor(ifreqs + ww);
ihi = ilow + 1;
ind = (ilow <= 0);
ilow(ind) = ilow(ind) + Nf;
ind = (ilow >= Nf+1);
ilow(ind) = ilow(ind) - Nf;
ind = (ihi >= Nf+1);
ihi(ind) = ihi(ind) - Nf;
for ir = 1:3
   for ic = 1:3
      icol = 3*(ir-1) + ic;
      SB(:,icol) = SB(:,icol)                                                    ...
                 + T0(ir,1)*((1-alf)*Spsi(ilow,1) + alf*Spsi(ihi,1))*Estar(ic,1) ...
                 + T0(ir,1)*((1-alf)*Spsi(ilow,2) + alf*Spsi(ihi,2))*Estar(ic,2) ...
                 + T0(ir,1)*((1-alf)*Spsi(ilow,3) + alf*Spsi(ihi,3))*Estar(ic,3) ...
                 + T0(ir,2)*((1-alf)*Spsi(ilow,4) + alf*Spsi(ihi,4))*Estar(ic,1) ...
                 + T0(ir,2)*((1-alf)*Spsi(ilow,5) + alf*Spsi(ihi,5))*Estar(ic,2) ...
                 + T0(ir,2)*((1-alf)*Spsi(ilow,6) + alf*Spsi(ihi,6))*Estar(ic,3) ...
                 + T0(ir,3)*((1-alf)*Spsi(ilow,7) + alf*Spsi(ihi,7))*Estar(ic,1) ...
                 + T0(ir,3)*((1-alf)*Spsi(ilow,8) + alf*Spsi(ihi,8))*Estar(ic,2) ...
                 + T0(ir,3)*((1-alf)*Spsi(ilow,9) + alf*Spsi(ihi,9))*Estar(ic,3);
   end
end

