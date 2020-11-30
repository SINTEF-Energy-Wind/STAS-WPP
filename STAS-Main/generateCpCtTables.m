clear;

%nm = 'Tjaereborg';
%lams = exp(linspace (log(2),log(45),50)).';
%bets = (pi/180)*[-5:45].';
%W = 2.5;

nm = 'DTU10MW';
lams = [0.1 [0.2:0.2:1] 1.2 1.5 [2:1:10] [12:2:20]].';
bets = (pi/180)*[[-4:2:18] 21 [25:5:90]].';
W = 1;

%lams = 10;
%bets = (pi/180)*90;

Nl = size (lams,1);
Nb = size (bets,1);

eval(['[s,a] = STASTurbine_' nm ' ();']);
load 'LTMnorms.txt'
length  = LTMnorms(1);
time    = LTMnorms(2);
mass    = LTMnorms(3);
current = LTMnorms(4);
power   = mass*(length^2)/(time^3);
voltage = sqrt(power);
velocity = length/time;
mass   = power*(time^3)/(length^2);
force  = mass*length/(time^2);
stress = force/(length^2);
ndens  = mass/(length^3);
nvisc  = mass/(length*time);
stiffness = force/length;
damping = force*time/length;
resistance = voltage/current;
inductance = voltage*time/current;
capacitance = current*time/voltage;
flux   = voltage*time;

dens = 1.225                           / ndens;
visc = 1.789e-5                        / nvisc;

psiFlag = 1;

Nel = a.Nb*a.Neb;

i2a = [1:2:2*Nel-1].';
i2b = [2:2:2*Nel].';
i3a = [1:3:3*Nel-2].';
i3b = [2:3:3*Nel-1].';
i3c = [3:3:3*Nel].';
i6a = [1:6:6*Nel-5].';
i6b = [2:6:6*Nel-4].';
i6c = [3:6:6*Nel-3].';
i6d = [4:6:6*Nel-2].';
i6e = [5:6:6*Nel-1].';
i6f = [6:6:6*Nel].';
i7a = [1:7:7*Nel-6].';
i7b = [2:7:7*Nel-5].';
i7c = [3:7:7*Nel-4].';
i7d = [4:7:7*Nel-3].';
i7e = [5:7:7*Nel-2].';
i7f = [6:7:7*Nel-1].';
i7g = [7:7:7*Nel].';

[idofs,idofm,inods,inodm,Ndof] = getDOFRefs (s);
[Tn_y,Th_d,Tb_h] = basicTransforms (s.nacelle.delta,s.driveshaft.phi);

Pin = assemblePin (s);

[qB0,Pn0_B,Ts0_B,TB0_g] = ...
         undeformedPosition (Pin,0,s.nacelle.delta,0,s.driveshaft.phi,[0;0;0],0, ...
                             idofs,idofm,inods,inodm);

t0 = 0;
Psi0 = 0;  % Must be zero here.

% Ignore structural pitch.
Tas = [0  0 -1; ...
      -1  0  0; ...
       0  1  0];

Trg = Tn_y;  % In this simplified case...

Tsp = zeros(3,3*Nel);
Tar = zeros(3,3*Nel);
Tpp0 = zeros(3,3*3);
Tp0g = zeros(3,3*3);
for ibod = 5:7

   jb3 = 3*(ibod-5);

   if (ibod == 5)
      idref = idofs(6);
      Neb   = s.blade(1).Nel;
      conns = s.blade(1).conn;
   elseif (ibod == 6)
      idref = idofs(7);
      Neb   = s.blade(2).Nel;
      conns = s.blade(2).conn;
   elseif (ibod == 7)
      idref = idofs(8);
      Neb   = s.blade(3).Nel;
      conns = s.blade(3).conn;
   end

   qB = qB0(idref+[1:6]);
   PB = Pn0_B(idref+[1:6]);

   Tpp0(:,jb3+[1:3]) = TFromTheta (qB(4:6));
   Tp0g(:,jb3+[1:3]) = TFromTheta (PB(4:6));

   ic3n = 3*Neb*(ibod-5);

   for iel = 1:Neb

      ic12 = 12*(iel-1);
      ic3  =  3*(iel-1);

      conn  = conns(:,iel);
      rdof  = idref + 6*(conn(1)-1);
      n1dof = idref + 6*(conn(2)-1);
      n2dof = idref + 6*(conn(3)-1);

      qn1 = qB0(n1dof+[1:6]);
      qn2 = qB0(n2dof+[1:6]);
      Pn1 = Pn0_B(n1dof+[1:6]);
      Pn2 = Pn0_B(n2dof+[1:6]);

      if (n1dof == rdof)

         qn1 = zeros(6,1);
         Pn1(1:3) = zeros(3,1);

         % Give the reference node the same undeformed orientation
         % as the second node.
         Pn1(4:6) = Pn0_B(n1dof+6+[4:6]);  

      end

      [xe,Tsp(:,ic3n+ic3+[1:3])] = elementCSFromNodes (qn1,qn2,Pn1,Pn2);

      Tar(:,ic3n+ic3+[1:3]) = (Trg.')*Tp0g(:,jb3+[1:3])*Tpp0(:,jb3+[1:3]) ...
                            * Tsp(:,ic3n+ic3+[1:3])*Tas;

   end

end

xer = zeros(3*Nel,1);
Lep = zeros(Nel,1);
xnr = zeros(3*(Nel+a.Nb),1);
xhg = globalPosition (qB0(idofs(4)+[1:6]),Pn0_B(idofs(4)+[1:6]), ...
                      qB0(idofm(6)-6+[1:6]),Pn0_B(idofm(6)-6+[1:6]));
Dia = 0;
for ib = 1:a.Nb

   ine = 3*a.Neb*(ib-1);
   inn = 3*(a.Neb+1)*(ib-1);

   qB = qB0(idofs(5+ib)+[1:6]);
   PB = Pn0_B(idofs(5+ib)+[1:6]);

   for ieb = 1:a.Neb
   
      i6 = 6*(ieb-1);
      i3 = 3*(ieb-1);

      qn2 = qB0(idofs(5+ib)+i6+[7:12]);
      Pn2 = Pn0_B(idofs(5+ib)+i6+[7:12]);
   
      if (ieb == 1)
         qn1 = zeros(6,1);
         Pn1 = zeros(6,1);
         Pn1(4:6) = Pn2(4:6);   
      else
         qn1 = qB0(idofs(5+ib)+i6+[1:6]);
         Pn1 = Pn0_B(idofs(5+ib)+i6+[1:6]);
      end

      xg1 = globalPosition (qB,PB,qn1,Pn1);
      xg2 = globalPosition (qB,PB,qn2,Pn2);
      xge = 0.5*(xg2 + xg1);

      xer(ine+i3+[1:3]) = (Trg.')*(xge - xhg);

      if (ieb == 1)
         xnr(inn+[1:3]) = (Trg.')*(xg1 - xhg);
      end
      xnr(inn+i3+[4:6]) = (Trg.')*(xg2 - xhg);

      Lep(a.Neb*(ib-1)+ieb) = sqrt((xnr(inn+i3+4) - xnr(inn+i3+1)).^2 ...
                            +      (xnr(inn+i3+5) - xnr(inn+i3+2)).^2);

   end

   val = sqrt(xnr(inn+3*(a.Neb-1)+4).^2 + xnr(inn+3*(a.Neb-1)+5).^2);
   Dia = Dia + val;

end

Dia = 2*Dia/3;
Area = (pi/4)*(Dia^2);

zr = zeros(2*Nel,1);
zr(i2a) = xer(i3a);
zr(i2b) = xer(i3b);

rs = zr(i2a);
r = [rs(1:Neb);rs(1:Neb);rs(1:Neb)];

%del = eps^0.5;
dell = 1e-2; % 1e-1;
delb = 1e-2;

ify = [2:6:6*Nel-4].';
ifz = [3:6:6*Nel-3].';
fid = fopen (['cpct_' nm '.txt'],'w');
for il = 1:Nl

%lam = 11.12014;
%bet3 = 0.06981317*[1;1;1];

   lam = lams(il);
   Vinf = (0.5*Dia*W/lam)*ones(Nel,1);

   for ib = 1:Nb

      bet3 = [bets(ib);bets(ib);bets(ib)];

      [Vir,Vmag,f,aoa,phi,Cl,Cd,dCl,dCd,Fld,Fas,Fap] =             ...
                         speedyBEM2 (a.Neb,Vinf,W,Psi0,bet3,a.xia, ...
                                     Dia,r,a.Lel,a.chord,a.xpc,    ...
                                     a.aoas,a.kfoils,a.foilwt,     ...
                                     dens,visc);
      Fz = -Fap(ify)*sin(bet3(1)) + Fap(ifz)*cos(bet3(1));
      Ft =  Fap(ify)*cos(bet3(1)) + Fap(ifz)*sin(bet3(1));
      Thrust = sum(Fz);
      Torq   = sum(r.*Ft);
      Pow    = W*Torq;
      Cpb = Pow/(0.5*dens*Area*(Vinf(1)^3));
      Ctb = Thrust/(0.5*dens*Area*(Vinf(1)^2));

      beth = bet3 + delb;
      [Vir,Vmag,f,aoa,phi,Cl,Cd,dCl,dCd,Fld,Fas,Fap] =             ...
                         speedyBEM2 (a.Neb,Vinf,W,Psi0,beth,a.xia, ...
                                     Dia,r,a.Lel,a.chord,a.xpc,    ...
                                     a.aoas,a.kfoils,a.foilwt,     ...
                                     dens,visc);
      Fz = -Fap(ify)*sin(beth(1)) + Fap(ifz)*cos(beth(1));
      Ft =  Fap(ify)*cos(beth(1)) + Fap(ifz)*sin(beth(1));
      Thrust = sum(Fz);
      Torq   = sum(r.*Ft);
      Pow    = W*Torq;
      Cpbh = Pow/(0.5*dens*Area*(Vinf(1)^3));
      Ctbh = Thrust/(0.5*dens*Area*(Vinf(1)^2));
      betl = bet3 - delb;
      [Vir,Vmag,f,aoa,phi,Cl,Cd,dCl,dCd,Fld,Fas,Fap] =             ...
                         speedyBEM2 (a.Neb,Vinf,W,Psi0,betl,a.xia, ...
                                     Dia,r,a.Lel,a.chord,a.xpc,    ...
                                     a.aoas,a.kfoils,a.foilwt,     ...
                                     dens,visc);
      Fz = -Fap(ify)*sin(betl(1)) + Fap(ifz)*cos(betl(1));
      Ft =  Fap(ify)*cos(betl(1)) + Fap(ifz)*sin(betl(1));
      Thrust = sum(Fz);
      Torq   = sum(r.*Ft);
      Pow    = W*Torq;
      Cpbl = Pow/(0.5*dens*Area*(Vinf(1)^3));
      Ctbl = Thrust/(0.5*dens*Area*(Vinf(1)^2));

      dCpdb = (Cpbh-Cpbl)/(2*delb);
      dCtdb = (Ctbh-Ctbl)/(2*delb);

      lamh = lam + dell;
      Vinfh = (0.5*Dia*W/lamh)*ones(Nel,1);
      [Vir,Vmag,f,aoa,phi,Cl,Cd,dCl,dCd,Fld,Fas,Fap] =              ...
                         speedyBEM2 (a.Neb,Vinfh,W,Psi0,bet3,a.xia, ...
                                     Dia,r,a.Lel,a.chord,a.xpc,     ...
                                     a.aoas,a.kfoils,a.foilwt,      ...
                                     dens,visc);
      Fz = -Fap(ify)*sin(bet3(1)) + Fap(ifz)*cos(bet3(1));
      Ft =  Fap(ify)*cos(bet3(1)) + Fap(ifz)*sin(bet3(1));
      Thrust = sum(Fz);
      Torq   = sum(r.*Ft);
      Pow    = W*Torq;
      Cplh = Pow/(0.5*dens*Area*(Vinfh(1)^3));
      Ctlh = Thrust/(0.5*dens*Area*(Vinfh(1)^2));
      laml = lam - dell;
      Vinfl = (0.5*Dia*W/laml)*ones(Nel,1);
      [Vir,Vmag,f,aoa,phi,Cl,Cd,dCl,dCd,Fld,Fas,Fap] =              ...
                         speedyBEM2 (a.Neb,Vinfl,W,Psi0,bet3,a.xia, ...
                                     Dia,r,a.Lel,a.chord,a.xpc,     ...
                                     a.aoas,a.kfoils,a.foilwt,      ...
                                     dens,visc);
      Fz = -Fap(ify)*sin(bet3(1)) + Fap(ifz)*cos(bet3(1));
      Ft =  Fap(ify)*cos(bet3(1)) + Fap(ifz)*sin(bet3(1));
      Thrust = sum(Fz);
      Torq   = sum(r.*Ft);
      Pow    = W*Torq;
      Cpll = Pow/(0.5*dens*Area*(Vinfl(1)^3));
      Ctll = Thrust/(0.5*dens*Area*(Vinfl(1)^2));

      dCpdl = (Cplh-Cpll)/(2*dell);
      dCtdl = (Ctlh-Ctll)/(2*dell);

      beth = bet3 + delb;
      betl = bet3 - delb;
      lamh = lam + dell;
      laml = lam - dell;
      Vinfh = (0.5*Dia*W/lamh)*ones(Nel,1);
      Vinfl = (0.5*Dia*W/laml)*ones(Nel,1);

      [Vir,Vmag,f,aoa,phi,Cl,Cd,dCl,dCd,Fld,Fas,Fap] =              ...
                         speedyBEM2 (a.Neb,Vinfh,W,Psi0,beth,a.xia, ...
                                     Dia,r,a.Lel,a.chord,a.xpc,     ...
                                     a.aoas,a.kfoils,a.foilwt,      ...
                                     dens,visc);
      Fz = -Fap(ify)*sin(beth(1)) + Fap(ifz)*cos(beth(1));
      Ft =  Fap(ify)*cos(beth(1)) + Fap(ifz)*sin(beth(1));
      Thrust = sum(Fz);
      Torq   = sum(r.*Ft);
      Pow    = W*Torq;
      Cphh = Pow/(0.5*dens*Area*(Vinfh(1)^3));
      Cthh = Thrust/(0.5*dens*Area*(Vinfh(1)^2));

      [Vir,Vmag,f,aoa,phi,Cl,Cd,dCl,dCd,Fld,Fas,Fap] =              ...
                         speedyBEM2 (a.Neb,Vinfh,W,Psi0,betl,a.xia, ...
                                     Dia,r,a.Lel,a.chord,a.xpc,     ...
                                     a.aoas,a.kfoils,a.foilwt,      ...
                                     dens,visc);
      Fz = -Fap(ify)*sin(betl(1)) + Fap(ifz)*cos(betl(1));
      Ft =  Fap(ify)*cos(betl(1)) + Fap(ifz)*sin(betl(1));
      Thrust = sum(Fz);
      Torq   = sum(r.*Ft);
      Pow    = W*Torq;
      Cphl = Pow/(0.5*dens*Area*(Vinfh(1)^3));
      Cthl = Thrust/(0.5*dens*Area*(Vinfh(1)^2)); 

      [Vir,Vmag,f,aoa,phi,Cl,Cd,dCl,dCd,Fld,Fas,Fap] =              ...
                         speedyBEM2 (a.Neb,Vinfl,W,Psi0,beth,a.xia, ...
                                     Dia,r,a.Lel,a.chord,a.xpc,     ...
                                     a.aoas,a.kfoils,a.foilwt,      ...
                                     dens,visc);
      Fz = -Fap(ify)*sin(beth(1)) + Fap(ifz)*cos(beth(1));
      Ft =  Fap(ify)*cos(beth(1)) + Fap(ifz)*sin(beth(1));
      Thrust = sum(Fz);
      Torq   = sum(r.*Ft);
      Pow    = W*Torq;
      Cplh = Pow/(0.5*dens*Area*(Vinfl(1)^3));
      Ctlh = Thrust/(0.5*dens*Area*(Vinfl(1)^2));

      [Vir,Vmag,f,aoa,phi,Cl,Cd,dCl,dCd,Fld,Fas,Fap] =              ...
                         speedyBEM2 (a.Neb,Vinfl,W,Psi0,betl,a.xia, ...
                                     Dia,r,a.Lel,a.chord,a.xpc,     ...
                                     a.aoas,a.kfoils,a.foilwt,      ...
                                     dens,visc);
      Fz = -Fap(ify)*sin(betl(1)) + Fap(ifz)*cos(betl(1));
      Ft =  Fap(ify)*cos(betl(1)) + Fap(ifz)*sin(betl(1));
      Thrust = sum(Fz);
      Torq   = sum(r.*Ft);
      Pow    = W*Torq;
      Cpll = Pow/(0.5*dens*Area*(Vinfl(1)^3));
      Ctll = Thrust/(0.5*dens*Area*(Vinfl(1)^2));

      d2Cp = (Cphh - Cphl - Cplh + Cpll)/(4*dell*delb);
      d2Ct = (Cthh - Cthl - Ctlh + Ctll)/(4*dell*delb);


      fprintf (fid,'%+5.6e %+5.6e %+5.6e %+5.6e %+5.6e %+5.6e %+5.6e %+5.6e %+5.6e %+5.6e\n', ...
               lam,bet3(1),Cpb,dCpdl,dCpdb,d2Cp,Ctb,dCtdl,dCtdb,d2Ct);

printf ('%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n', ...
        lam,bet3(1),Cpb,dCpdl,dCpdb,d2Cp,Ctb,dCtdl,dCtdb,d2Ct);

%{

% Numerically unstable derivatives at certain conditions, likely a result of fsolve () in speedyBEM2.

      betc = bet3 + i*del;
      [Vir,Vmag,f,aoa,phi,Cl,Cd,dCl,dCd,Fld,Fas,Fap] =             ...
                         speedyBEM2 (a.Neb,Vinf,W,Psi0,betc,a.xia, ...
                                     Dia,r,a.Lel,a.chord,a.xpc,    ...
                                     a.aoas,a.kfoils,a.foilwt,     ...
                                     dens,visc);
      Fz = -Fap(ify)*sin(betc(1)) + Fap(ifz)*cos(betc(1));
      Ft =  Fap(ify)*cos(betc(1)) + Fap(ifz)*sin(betc(1));
      Thrust = sum(Fz);
      Torq   = sum(r.*Ft);
      Pow    = W*Torq;
      Cpb = Pow/(0.5*dens*Area*(Vinf(1)^3));
      Ctb = Thrust/(0.5*dens*Area*(Vinf(1)^2));

      lamc = lam + i*del;
      Vinfc = (0.5*Dia*W/lamc)*ones(Nel,1);
      [Vir,Vmag,f,aoa,phi,Cl,Cd,dCl,dCd,Fld,Fas,Fap] =              ...
                         speedyBEM2 (a.Neb,Vinfc,W,Psi0,bet3,a.xia, ...
                                     Dia,r,a.Lel,a.chord,a.xpc,     ...
                                     a.aoas,a.kfoils,a.foilwt,      ...
                                     dens,visc);
      Fz = -Fap(ify)*sin(bet3(1)) + Fap(ifz)*cos(bet3(1));
      Ft =  Fap(ify)*cos(bet3(1)) + Fap(ifz)*sin(bet3(1));
      Thrust = sum(Fz);
      Torq   = sum(r.*Ft);
      Pow    = W*Torq;
      Cpl = Pow/(0.5*dens*Area*(Vinfc(1)^3));
      Ctl = Thrust/(0.5*dens*Area*(Vinfc(1)^2));

      fprintf (fid,'%+5.6e %+5.6e %+5.6e %+5.6e %+5.6e %+5.6e %+5.6e %+5.6e\n', ...
               lam,bet3(1),real(Cpb),imag(Cpl)/del,imag(Cpb)/del,               ...
                           real(Ctb),imag(Ctl)/del,imag(Ctb)/del);

%printf ('%12.4f %12.4f %12.4f %12.4f\n',lam,bet3(1),real(Cpb),real(Ctb));
printf ('%12.4f %12.4f %12.4f %12.4f %12.4f\n',lam,bet3(1),real(Cpb),imag(Cpl)/del,imag(Cpb)/del);

%}


fflush(stdout);

%return

   end


end

fclose('all');


