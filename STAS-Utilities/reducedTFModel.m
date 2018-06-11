function [Ar,Bur,Cr,Dur] = reducedTFModel (mflag,A,B,C,D,fs,acc)


% CAUTION, NOT ADAPTED FOR COMPLEX STEP DERIVATIVES.


%
% Generate a reduced model from a full model.  Two methods are
% implemented.  One is modal truncation.  For a given set of
% control-to-sensor transfer functions, a set of modes is
% determined that reproduce these transfer functions to a given
% relative accuracy.  The second method is a Hankel-norm method
% based on the control toolbox.
%
% Version:        Changes:
% --------        -------------
% 08.05.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 08.05.2017      
%
% Inputs:
% -------
% mflag           : 1 for modal truncation, 2 for bal. realization.
% A               : Full system matrix.
% B               : Control excitation of state equations.
% C               : Sensor output dependence on states.
% D               : Sensor output dependence directly on controls.
% fs              : Modes: Frequencies at which the TFs shall match to
%                   the specified accuracy.
% acc             : Modes: relative accuracy.  Hankel: number of
%                   states in the system.  
%
% Outputs:
% --------
%

Nx = size(A,1);
Nu = size(B,2);
Ny = size(C,1);
Nf = size(fs,1);

if (mflag == 1)

   % Modes.  slap contains the eigenvalues in order of increasing
   % frequency (imaginary part), shp the eigenvectors in each
   % column, and ifrq the column indices of shp that correspond
   % to the sorted eigenvalues of slap.
   [slap,shp,ifrq] = eigVal (A);
   Nmod = size(shp,2);
   shn = zeros(size(shp));
   for imod = 1:Nmod
      shmag = sqrt(shp(:,ifrq(imod))'*shp(:,ifrq(imod)));
      shn(:,imod) = shp(:,ifrq(imod)) ... % Sort and normalize.
                  / shmag;
   end
   shp = shn;
   ii  = 1:Nx;
   jj  = ii;
   Lam = sparse(ii,jj,slap,Nx,Nx);      % Orthogonalized A matrix.

   % Identify the pairs of complex conjugate modes.
   Npairs = sum(abs(imag(slap)) > 0)/2; % Number of complex-conjugate modes.
   pairs  = [[1:Npairs].' [Nx:-1:Nx-Npairs+1].'];
   Nsp    = Nx - 2*Npairs;              % Number of single-pole modes.
   Nmod   = Nsp + Npairs;

   ss  = ones(Nx,1);
   INx = sparse(ii,jj,ss,Nx,Nx);        % Identity matrix.

   res  = zeros(Nmod*Nf*Ny,Nu);
   ishp = balanceDiv (shp,eye(size(shp)));
   ishB = ishp*B;
   Csh  = C*shp;
   mmet = zeros(Nmod,1);
   for jf = 1:Nf

      w = 2*pi*fs(jf);

      % The full transfer functions at the selected frequencies.
      IL = i*w*INx - Lam;
      dydu = Csh*(IL\ishB); % + D;  Don't include D in the evaluation.

      % The modal projections, one value for each mode, at each
      % frequency, for each control-to-sensor TF.  Also compute
      % metrics of influence.
      for ipair = 1:Npairs
         rref = Nmod*Ny*(jf-1) + Ny*(ipair-1);

         for ip = 1:2
            imod = pairs(ipair,ip);
            TFmod = Csh(:,imod)*ishB(imod,:)/(i*w - slap(imod));
% A problem with using real() is that then although the resulting
% projection of the TF accurately predicts the actual TF, there
% may still be a significant imaginary part that makes the 
% reduced TF inaccurate, in both magnitude and phase.
%            res(rref+[1:Ny],:) = res(rref+[1:Ny],:)      ...
%                               + real(conj(TFmod).*dydu) ...
%                              ./ (abs(dydu).^2);
            res(rref+[1:Ny],:) = res(rref+[1:Ny],:)      ...
                               + conj(TFmod).*dydu       ...
                              ./ (abs(dydu).^2);
         end

% I can sort via mmet according to real() here if desired, such
% that the modes are selected in the order of greatest contribution
% to the projection in the direction of the true TF, as opposed to
% the raw magnitude.
%         mmet(ipair) = mmet(ipair) ...
%                     + sum(sum(abs(real(res(rref+[1:Ny],:))),2));
         mmet(ipair) = mmet(ipair) ...
                     + sum(sum(abs(res(rref+[1:Ny],:)),2));
      end

      for isp = 1:Nsp
         rref = Nmod*Ny*(jf-1) + Ny*Npairs + Ny*(isp-1);
         imod = Npairs + isp;
         TFmod = Csh(:,imod)*ishB(imod,:)/(i*w - slap(imod));
         res(rref+[1:Ny],:) = conj(TFmod).*dydu ...
                           ./ (abs(dydu).^2);
%         mmet(imod) = mmet(imod) ...
%                    + sum(sum(abs(real(res(rref+[1:Ny],:))),2));
         mmet(imod) = mmet(imod) ...
                    + sum(sum(abs(res(rref+[1:Ny],:)),2));
      end

   end

   % Sort according to the metrics of influence.
   [srt,sind] = sort(mmet,'descend');  

   % Determine which modes to retain.
   resc = zeros(Ny*Nf,Nu);
   imod = 0;
   magres = 0;
   phres = 0;
   while ((imod == 0) || (((magres > acc) || (phres > acc)) && (imod < Nmod)))
      imod++;
      for jf = 1:Nf         
         rref = Nmod*Ny*(jf-1) + Ny*(sind(imod)-1);
         resc(Ny*(jf-1)+[1:Ny],:) = resc(Ny*(jf-1)+[1:Ny],:) ...
                                  + res(rref+[1:Ny],:);
      end
[imod sind(imod)]
resc
      magres = max(max(max(abs(resc)-1,1-abs(resc))));
      phres  = max(max(abs(atan2(imag(resc),real(resc))/pi)));

[magres phres]

max(magres,phres)

   end
   Nret = imod;                          % Keep the first Nret entries
                                         % in sind.

   % Build the reduced matrices.  This involves transforming complex
   % conjugate modes into real.  In order to do this, we must work
   % with modal pairs.  The transform applied to the Lambda matrix
   % for a complex conjugate pair of modes is Y = [1 1;-i i].  This
   % premultiplies the Lambda and B matrices, and its inverse
   % postmultiplies Lambda.
   % 
   % The list sind includes the first mode of each complex pair, and
   % the real-root modes.  If therefore the sorted index is less 
   % than Npairs, this means that the pair of conjugates is to be
   % included in Lambda.
   Ncm = sum(sind(1:Nret)<=Npairs);   % Number of retained complex modes.
   Nre = Nret - Ncm;                  % Number of retained real modes.
   Nst = 2*Ncm + Nre;                 % Number of states retained.
   ii = [1:Nst];
   jj = [1:Nst];
   ss = ones(Nst,1);
   YY = sparse(ii,jj,ss,Nst,Nst);     % Initialize with the identity matrix.
   mat = [1 1;-i i];

   si = zeros(Nst,1);
   idof = 0;
   for imod = 1:Nret

      if (sind(imod) <= Npairs)
         % For each complex mode, replace the values in YY with the
         % desired transform.  The modes will be ordered such that the
         % complex pairs appear next to each other in the state vector.
         YY(idof+[1:2],idof+[1:2]) = mat;

         % Also collect the state indices, while I'm at it.  This is
         % like sind, but includes all the pairs of complex modes, not
         % just the first.
         si(idof+[1:2]) = [sind(imod) Nx-sind(imod)+1].';

         idof = idof + 2;
      else

         % Real mode.
         si(idof+1) = sind(imod);

         idof++;

      end

   end
   iYY = inv(YY);
[[1:size(si,1)].' si real(slap(si))/(2*pi) imag(slap(si))/(2*pi)]
   % Define the baseline Lambda matrix.
   ii = [1:Nst];
   jj = ii;
   ss = slap(si);
   Lambda = sparse(ii,jj,ss,Nst,Nst);

%'Get rid of this scaling when testLQG9 is fixed'
%scl = diag(1./shp(177,si))
%iscl = inv(scl)

   % Transform to the Ar matrix.
   Ar = real(YY*Lambda*iYY);
%   Ar = real(YY*iscl*Lambda*scl*iYY);

   % Transform and reorganize the Br and Cr matrices.  real() is to
   % get rid of minute numerical noise in the complex portion of the
   % arrays.
   Bur = real(YY*ishB(si,:));
   Cr = real(Csh(:,si)*iYY);
%   Bur = real(YY*iscl*ishB(si,:));
%   Cr = real(Csh(:,si)*scl*iYY);

   Dur = D;

elseif (mflag == 2)

   GG = ss (A,B,C,D);
   [Gr,info] = hnamodred (GG,acc);
   [Ar,Bur,Cr,Dur] = ssdata (Gr);

end


