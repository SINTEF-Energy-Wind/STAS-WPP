clear;

nm = 'DTU10MW';
eval(['[s,a] = STASTurbine_'  nm ' ();']);

Nmud = s.foundation.Nmud;
Nwater = s.foundation.Nwater;
Lel = s.foundation.Lel;
Diaf = s.foundation.D;

depth = sum(Lel(Nmud+1:Nmud+Nwater));
Nnod = s.foundation.Nnod - Nmud;
znod = zeros(Nnod,1);  
znod(1) = -depth;
for inod = 2:Nnod
   znod(inod) = znod(inod-1) + Lel(Nmud+inod-1);
end
Dnod = zeros(Nnod,1);
Dnod(1) = Diaf(Nmud+1);
Dnod(2:Nnod) = Diaf(Nmud+1:Nmud+Nnod-1);

%{

% The coarse Hs-Tp table from Memo AN 15.12.19.
Hs = [1;1;1;1;1;1;1;1;2;2;2;2;2;2;2;3;3;3;3;3; ...
      4;4;4;4;4;5;5;5;5;6;6;6;7;7;7;8];
Tp = [2;4;6;8;10;12;14;16;4;6;8;10;12;14;16;6;8;10;12;14; ...
      6;8;10;12;14;8;10;12;14;10;12;14;10;12;14;14];

%}
%{

% A refined Hs-Tp table.
Hs = [0.5*ones(16,1);1.0*ones(18,1);1.5*ones(17,1);2.0*ones(16,1); ...
      2.5*ones(14,1);3.0*ones(13,1);3.5*ones(12,1);4.0*ones(11,1); ...
      4.5*ones(10,1);5.0*ones(10,1);5.5*ones(9,1);6.0*ones(7,1);   ...
      6.5*ones(7,1);7.0*ones(6,1);7.5*ones(5,1);8.0*ones(4,1);     ...
      8.5*ones(4,1);9.0*ones(3,1)];
Tp = 1 + [[1:16]';[1:18]';[2:18]';[2:17]';   ...
          [3:16]';[4:16]';[4:15]';[5:15]';   ...
          [6:15]';[6:15]';[7:15]';[8:14]';   ...
          [8:14]';[9:14]';[10:14]';[11:14]'; ...
          [11:14]';[12:14]'];

%}




% A single analysis...
Hs = 4.0;
Tp = 8.0;
Nf = 2^12;
df = 0.001;
rhow = 1025;


Nwave = size(Hs,1);

% Define elements for the surface pressure analysis, with the z
% coordinate relative to the still water line.  The coordinates
% specify the element centroids.  Theta is measured relative to
% the direction of wave propagation.
z = [[-30:2:-16] [-14.5:1:-5.5] [-4.75:0.5:-0.25]]';
theta = (pi/180)*[0:10:180]';

Nz = size(z,1);
Vc = zeros(2,Nz);
Vc(1,:) = 0;
Vc(2,:) = 0;

T = 1/df;
dt = 1/(Nf*df);
f = [df:df:(Nf/2)*df]';  % For sampling the one-sided wave spectrum.
t = [0:dt:(Nf-1)*dt]';

Niter = 20;
for iw = 1:Nwave

printf('Case %6d of %6d\n',iw,Nwave);
fflush(stdout);

   % Take the average over a fairly large number of iterations in order
   % to smooth the spectrum.  Random numbers are employed in these
   % calculations.
   for iter = 1:Niter

printf('Iteration %6d of %6d\n',iter,Niter);
fflush(stdout);

      F = waveForce (Hs(iw),Tp(iw),Vc,depth,rhow,f,theta,z,t,znod,Dnod);

% Fsum can be used for accumulating the forces to, for instance, the
% waterline, for simplified analysis.
%      Fsum = zeros(Nf,1);
%      for inod = 1:Nnod
%         iref = Nf*(inod-1);
%         Fsum = Fsum + F(iref+1:iref+Nf);
%      end

      [Qij,Sij] = forceSpectra (F,Nf,dt,Nnod);
%      [Qij,Sij] = forceSpectra (Fsum,Nf,dt,1);

      if (iter == 1)
         Qavg = Qij;
         Savg = Sij;
      else
         Qavg = (Qavg + Qij/(iter-1))*(iter-1)/iter;
         Savg = (Savg + Sij/(iter-1))*(iter-1)/iter;
      end

   end

%{

   % Additional smoothing may be employed, if desired, by uncommenting
   % this for loop block.
   for inod1 = 1:Nnod
      for inod2 = 1:Nnod
         ind = Nnod*(inod1-1) + inod2;

         % Do not smooth the zero frequency (mean).  Otherwise,
         % apply a five-point moving average.
         Smoov = zeros(Nf,1);
         Smoov(1) = Savg(1,ind);
         Smoov(2) = 0.2*(3*Savg(2,ind) + Savg(3,ind) + Savg(4,ind));
         Smoov(3) = 0.2*(2*Savg(2,ind) + Savg(3,ind) ...
                  +      Savg(4,ind) + Savg(5,ind));
         Smoov(Nf-1) = 0.2*(Savg(Nf-3,ind) + Savg(Nf-2,ind) ...
                     +      Savg(Nf-1,ind) + 2*Savg(Nf,ind));
         Smoov(Nf) = 0.2*(Savg(Nf-2,ind) + Savg(Nf-1,ind) ...
                   +      3*Savg(Nf,ind));
         Smoov(4:Nf-2) = 0.2*(Savg(2:Nf-4,ind) ...
                       +      Savg(3:Nf-3,ind) ...
                       +      Savg(4:Nf-2,ind) ...
                       +      Savg(5:Nf-1,ind) ...
                       +      Savg(6:Nf,ind));

         Savg(:,ind) = Smoov;

      end
   end

%}

   save('-binary',['Savg_Hs' int2str(10*Hs(iw)) '_TP' ...
                   int2str(Tp(iw)) '.bin'],'Savg');

end

%{

% Some possible outputs ...

inod1 = 6;
iref1 = Nf*(inod1-1);
inod2 = 11;
iref2 = Nf*(inod2-1);
fid3 = fopen('F.txt','w');
for it = 1:Nf
%   fprintf(fid3,'%+5.4e %+5.4e %+5.4e\n', ...
%           t(it),F(iref1+it),F(iref2+it));
   fprintf(fid3,'%+5.4e %+5.4e\n', ...
           t(it),Fsum(it));
end

fid = fopen('Q.txt','w');
inod1 = 1; % 11;
inod2 = 1; % 11;
ind = Nnod*(inod1-1) + inod2;
for it = 1:Nf
   fprintf(fid,'%+5.4e %+5.4e\n',t(it),real(Qavg(it,ind)));
end

fid2 = fopen('S.txt','w');
ind = 1; % Nnod*(inod1-1) + inod2;
for ifreq = 1:Nf
   fprintf(fid2,'%+5.4e %+5.4e\n',df*(ifreq-1),real(Savg(ifreq,ind)));
end

fclose('all');

%}


