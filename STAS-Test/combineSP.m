function nSP = combineSP (sbin,nS,sP,nalf)
%
% Combine cycle counts of a stochastic signal s_S with a periodic
% signal s_P.
%
% Version:        Changes:
% --------        -------------
% 17.09.2020      Original code.
%
% Version:        Verification:
% --------        -------------
% 17.09.2020      
%
% Inputs:
% -------
% sbin            : [smin,ds,smax]
% nS              : Cycle counts for each bin.
% sP              : Amplitude of periodic signal.
% nalf            : Number of relative phases to consider.
%
% Outputs:
% --------
% 

n = size(nS,1);
nSP = zeros(n,1);

for ia = 1:nalf
%   alf = 2*pi*(ia-1)/nalf;
   alf = pi*(ia-1)/(nalf-1);  % Take advantage of symmetry.
   sb = sbin(1) + sbin(2)*([0:n-1]);
   AA = sb + sP*cos(alf);
   BB = -sP*sin(alf);
   th = atan2(BB,AA);
   s = abs(AA.*cos(th) + BB.*sin(th));
   [ic,wc] = nearestCells (s,sbin);

   % Leverage sparse to sum repeated indices.
   vals = nS.*(wc(1,:).')/nalf;
   vec = sparse(ic(1,:),1,vals,n,1);
   nSP = nSP + vec;
   vals = nS.*(wc(2,:).')/nalf;
   vec = sparse(ic(2,:),1,vals,n,1);
   nSP = nSP + vec;

% Incorrect because this overwrites repeated indices, rather than summing.
%   nSP(ic(1,:)) = nSP(ic(1,:)) + nS.*(wc(1,:).')/nalf;
%   nSP(ic(2,:)) = nSP(ic(2,:)) + nS.*(wc(2,:).')/nalf;

%'----------------------'
%alf
%full([ic;wc])
%full([sb.' s.' nS wc(1,:).' wc(2,:).' nS.*(wc(1,:).')/nalf nS.*(wc(2,:).')/nalf])

end