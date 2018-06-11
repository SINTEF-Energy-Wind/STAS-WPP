function [shape,freq,mdamp] = bodyModes (s,M,K,ret,slv)
%
% Version:        Changes:
% --------        -------------
% 04.12.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 04.12.2017      Visual check of all modes.  Experience from similar
%                 legacy code.
%

%'bodyModes'

Ndof = size(M,1);

[idofs,idofm,inods,inodm,Nd0] = getDOFRefs (s);

TotDOFs = 6;
TotMods = 12;  % 6 foundation ref DOFs, 6 joint DOFs.
for ibod = 1:7
   if (ibod == 1)
      Nnod = s.foundation.Nnod;
      Nmod = s.foundation.Nmod;
   elseif (ibod == 2)
      Nnod = s.tower.Nnod;
      Nmod = s.tower.Nmod;
   elseif (ibod == 3)
      Nnod = s.nacelle.Nnod;
      Nmod = s.nacelle.Nmod;
   elseif (ibod == 4)
      Nnod = s.driveshaft.Nnod;
      Nmod = s.driveshaft.Nmod;
   elseif (ibod == 5)
      Nnod = s.blade(1).Nnod;
      Nmod = s.blade(1).Nmod;
   elseif (ibod == 6)
      Nnod = s.blade(2).Nnod;
      Nmod = s.blade(2).Nmod;
   elseif (ibod == 7)
      Nnod = s.blade(3).Nnod;
      Nmod = s.blade(3).Nmod;
   end
   TotMods = TotMods + Nmod;
   TotDOFs = TotDOFs + 6*Nnod*Nmod;
end

shape = spalloc(Ndof,TotMods,TotDOFs);
freq  = zeros (TotMods,1);
mdamp = zeros (TotMods,1);

shape(1:6,1:6) = eye(6);  % Foundation ref nodes.  Zero frequency.
imod = 6;
ind = [1:Ndof].';
for ibod = 1:7

   if (ibod == 1)
      Nmod = s.foundation.Nmod;
      DOF1 = ind(ret==idofs(1)+7);
      DOF2 = ind(ret==idofs(2));
   elseif (ibod == 2)
      Nmod = s.tower.Nmod;
      DOF1 = ind(ret==idofs(2)+7);
      DOF2 = ind(ret==idofs(3));
   elseif (ibod == 3)
      Nmod = s.nacelle.Nmod;
      DOF1 = ind(ret==idofs(3)+7);
      DOF2 = ind(ret==idofs(4));
   elseif (ibod == 4)
      Nmod = s.driveshaft.Nmod;
      DOF1 = ind(ret==idofs(4)+7);
      DOF2 = ind(ret==idofs(6));
   elseif (ibod == 5)
      Nmod = s.blade(1).Nmod;
      DOF1 = ind(ret==idofs(6)+7);
      DOF2 = ind(ret==idofs(7));
   elseif (ibod == 6)
      Nmod = s.blade(2).Nmod;
      DOF1 = ind(ret==idofs(7)+7);
      DOF2 = ind(ret==idofs(8));
   elseif (ibod == 7)
      Nmod = s.blade(3).Nmod;
      DOF1 = ind(ret==idofs(8)+7);
      DOF2 = ind(ret==idofs(8)+(idofs(7)-idofs(6)));
   end

   % Force symmetry in the mass and stiffness matrices used to calculate
   % the body modes.  Small antisymmetric gyroscopic terms have been
   % observed to produce complex mode shapes for some selected modes,
   % which then destroys the real state space calculation.  Note that
   % I am not using these modes in any sense that requires that they
   % precisely diagonalize M and K.  It's OK if they are approximate.
   %
   % Also, it is vital that the order of the modes is the same for the
   % three blades.  Enforce this by requiring the blade 1 mode shapes
   % to be used for each blade.
   dofs = [DOF1:DOF2].';
   if (ibod <= 5)
      % Compute body modes.  ibod = 6 or 7, use modes from ibod = 5.
      Msym = 0.5*(M(dofs,dofs) + M(dofs,dofs).');
      Ksym = 0.5*(K(dofs,dofs) + K(dofs,dofs).');
      [shp,frq] = modes (Nmod,Msym,Ksym);
   end
   shape(dofs,imod+[1:Nmod]) = shp;
   freq(imod+[1:Nmod]) = frq;

   if (s.zdamp > 0)
      Mmod = (shp.')*Msym*shp;
      Kmod = (shp.')*Ksym*shp;
      mdamp(imod+[1:Nmod]) = 2*s.zdamp ...
                           * sqrt(absc(diag(Mmod).*diag(Kmod)));
   end

%{
%if (ibod == 1)
ibod
   diag(M(dofs,dofs))
   diag(K(dofs,dofs))
   %shp
   frq
nd = size(dofs,1);
fidm = fopen('M.txt','w');
fidk = fopen('K.txt','w');
for ir = 1:nd
   fprintf(fidm,'%+5.6e',M(dofs(ir),dofs(1)));
   fprintf(fidk,'%+5.6e',K(dofs(ir),dofs(1)));
   for ic = 2:nd
      fprintf(fidm,' %+5.6e',M(dofs(ir),dofs(ic)));
      fprintf(fidk,' %+5.6e',K(dofs(ir),dofs(ic)));
   end
   fprintf(fidm,'\n');
   fprintf(fidk,'\n');
end
fclose('all');
%end
%}

%{
'==============================='
ibod
for jj = 1:Nmod
   shf = zeros(size(ret)+size(slv),1);
   shf(ret) = shape(:,imod+jj);
   freq(imod+jj)
   printVec(shf)
   ' '
end
%}

   imod = imod + Nmod;

end

shape(DOF2+[1:6],imod+[1:6]) = eye(6);  % Joint DOFs.  Zero frequency.



