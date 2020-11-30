function [Ar,Br,Cr,Dr] = modelReduction (A,B,C,D,freqs,tol)
%
% Creates a reduced-order model that satisfies all transfer functions
% to within a given tolerance, at a set of specified frequencies.
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
% A, B, C, D      : State matrices.
% freqs           : Frequencies at which the TF match is evaluated.
% tol             : Absolute tolerance to which the transfer functions
%                   must be satisfied. 
%
% Outputs:
% --------
% Ar, Br, Cr, Dr  : State matrices.

Nx = size(A,1);
Nu = size(B,2);
Ny = size(C,1);
Nf = size(freqs,1);
Ntf = Ny*Nu

% Evaluate the full TFs at the specified frequencies.
Hf = zeros(Ntf,Nf);
for ifreq = 1:Nf
   w = 2*pi*freqs(ifreq);
   irow = Ntf*(ifreq-1);
   TFmat = C*((i*w*speye(Nx) - A)\B) + D;
   Hf(:,ifreq) = reshape (TFmat,Ntf,1);
end

G = ss (A,B,C,D);

Nxr = 2;
conv = 0;
while (conv == 0)

   Nxr = Nxr + 1;

   [Gr,outp] = spamodred (G,Nxr);
   Ar = Gr.a;
   Br = Gr.b;
   Cr = Gr.c;
   Dr = Gr.d;

   H = zeros(Ntf,Nf);
   for ifreq = 1:Nf
      w = 2*pi*freqs(ifreq);
      irow = Ntf*(ifreq-1);
      TFmat = Cr*((i*w*speye(Nxr) - Ar)\Br) + Dr;
      H(:,ifreq) = reshape (TFmat,Ntf,1);
   end

   err = max(max(abs(H - Hf)));

printf ('%5d %+5.6e\n',Nxr,err);
fflush(stdout);

   if ((err <= tol) || (Nxr == Nx))
      conv = 1;
   end
      
   end

end

