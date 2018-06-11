function [Qij,Sij] = forceSpectra (F,Nf,dt,Nnod)
%
% The wave forces are computed at the given nodal locations.
% The cross-spectral matrix between the wave forces at the
% nodes is computed here.
%
% Version:        Changes:
% --------        -------------
% 25.03.2015      Original code.
%
% Version:        Verification:
% --------        -------------
% 25.03.2015      Memo AN 15.12.19
%
% Inputs:
% -------
% F               : Time history of wave-direction forces on each
%                   node.  Vector of length Ntime*Nnod.
% Nf              : Number of timesteps/number of frequencies/
%                   number of time offsets.
% dt              : Timestep, for scaling the FFT.
% Nnod            : Number of nodes for cross-correlations.
%
% Outputs:
% --------
% Qij,Sij         : Correlations and spectra.

Sij = sparse(Nf,Nnod*Nnod);
Qij = sparse(Nf,Nnod*Nnod);
FF  = sparse(Nf,Nnod);

for ii = 1:Nnod
   indi = Nf*(ii-1);
   Fi = zeros(Nf,1);
   Fi(1:Nf/2) = F(indi+Nf/2+1:indi+Nf);
   Fi(Nf/2+1:Nf) = F(indi+1:indi+Nf/2);
   FF(:,ii) = fft(Fi);
end

for ii = 1:Nnod
   ind = Nnod*(ii-1);

   for jj = 1:Nnod

      % Compute the correlation between the two force vectors
      % using FFTs.  (Press et al. pp 602 and 648.)
      prod = FF(:,jj).*conj(FF(:,ii))/Nf;
      Sij(:,ind+jj) = dt*prod;
      Qij(:,ind+jj) = ifft(prod);

   end

end
