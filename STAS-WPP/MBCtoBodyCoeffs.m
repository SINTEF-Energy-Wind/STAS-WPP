function [cB,rB] = MBCtoBodyCoeffs (opt,cpin,rpsi,azi,W,df)
%
% Convert an input set of Fourier coefficients in multi-blade or
% nacelle coordinates to an output set in rotating coordinates.
% Options are provided for transforming from multi-blade coordinates
% to each blade, and from nacelle to driveshaft coordinates.
%
% The function also transforms periodic Fourier coefficients and
% provides the mean values in body coordinates.
%
% If an MBC transform is being performed then the indices represent
% the blade number, and if a driveshaft transform is being performed
% then the indices represent 1=x, 2=y, 3=z.
%
% Version:        Changes:
% --------        -------------
% 08.09.2020      Original code.
%
% Version:        Verification:
% --------        -------------
% 08.09.2020      Verified using an independent calculation based
%                 on a time-domain transform of the periodic signal.
%
% Inputs:
% -------
% opt             : = 1 for MBC to body,
%                   = 2 for driveshaft.
% cpin            : Fourier coefficients representing a periodic signal.
%                   rper(p) = sum_k c_k exp(-i 2 pi (k p / N)).
% rpsi            : mean in MBC coordinates.
% azi             : azimuth angle at t = 0 (usually 0).
% W               : Rotor speed, rad/s.
% df              : Frequency bin, Hz.
%
% Outputs:
% --------
% cB              : periodic coefficients in body coordinates.
% rB              : the mean in body coordinates.

Nf = size(cpin,1);
ww = W/(2*pi*df);
cpsi = cpin;

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

% Mean value.
ind = round([Nf-ww+1;ww+1]);
jj = (ind > Nf);
ind(jj) = ind(jj) - Nf + 1;

rB = E0*rpsi + E*cpsi(ind(1),:).' + Estar*cpsi(ind(2),:).';

% Coefficients.  Initialize with the E0 terms.
cB = [E0(1,1)*cpsi(:,1)+E0(1,2)*cpsi(:,2)+E0(1,3)*cpsi(:,3), ...
      E0(2,1)*cpsi(:,1)+E0(2,2)*cpsi(:,2)+E0(2,3)*cpsi(:,3), ...
      E0(3,1)*cpsi(:,1)+E0(3,2)*cpsi(:,2)+E0(3,3)*cpsi(:,3)];

ifreqs = [1:Nf].';

% Augment the c0 term with rpsi, which applies from here on, but
% not before.  This allows the calculations below to be done in
% one consistent step.
cpsi(1,:) = rpsi.';

% The E c_(k-w) term.
indw = round(ifreqs - ww);
ind = (indw <= 0);
indw(ind) = indw(ind) + Nf;
cB = cB + [E(1,1)*cpsi(indw,1) + E(1,2)*cpsi(indw,2) + E(1,3)*cpsi(indw,3), ...
           E(2,1)*cpsi(indw,1) + E(2,2)*cpsi(indw,2) + E(2,3)*cpsi(indw,3), ...
           E(3,1)*cpsi(indw,1) + E(3,2)*cpsi(indw,2) + E(3,3)*cpsi(indw,3)];

% The E* c_(k+w) term.
indw = round(ifreqs + ww);
ind = (indw >= Nf+1);
indw(ind) = indw(ind) - Nf;
cB = cB + [Estar(1,1)*cpsi(indw,1) + Estar(1,2)*cpsi(indw,2) + Estar(1,3)*cpsi(indw,3), ...
           Estar(2,1)*cpsi(indw,1) + Estar(2,2)*cpsi(indw,2) + Estar(2,3)*cpsi(indw,3), ...
           Estar(3,1)*cpsi(indw,1) + Estar(3,2)*cpsi(indw,2) + Estar(3,3)*cpsi(indw,3)];

% The imperfect binning and rounding of the rotor frequency results in
% adjacent frequency bins filled with values that should in reality
% overlap and cancel.  Fix this by summing to the most-central frequency.
cB(1,:) = 0.;  % By definition, see derivations.
cB(Nf,:) = 0.; % Also zero the adjacent bins.
cB(2,:) = 0.;

jf = ww;
while (jf <= Nf/2-1)

   jj = round(jf) + 1;
   cB(jj,:) = cB(jj-1,:) + cB(jj,:) + cB(jj+1,:);
   cB(jj-1,:) = 0.;
   cB(jj+1,:) = 0.;

   jf = jf + ww;

end % while

jf = Nf - ww + 1;
while (jf >= Nf/2)

   jj = round(jf);
   cB(jj,:) = cB(jj-1,:) + cB(jj,:) + cB(jj+1,:);
   cB(jj-1,:) = 0.;
   cB(jj+1,:) = 0.;

   jf = jf - ww;

end % while
