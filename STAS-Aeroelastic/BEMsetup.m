function [Tas,ch,Lel,foilwt,aoaz,aoast,            ...
          xas,yas,iq] = BEMsetup (s,a)
%
% Define some setup parameters for the calls to BEMNL or BEMlin.
%
% Version:        Changes:
% --------        -------------
% 03.03.2018      Original code.
%
% Version:        Verification:
% --------        -------------
% 03.03.2018      
%
% Inputs:
% -------
% s,a             : Structures containing structural and aerodynamic
%                   definitions.
% P               : Undeformed nodal positions.
% aoas,kfoils     : Splined airfoil coefficient tables.
%
% Outputs:
% --------
% Tas             : Transform from airfoil to section coordinates.
%                   3-by-3*Nel.
% TB0g            : Transform from undeformed body to global coordinates.
%                   3-by-3*Nel.
% ch,Lel          : Chord and lengths of the elements.  Length is measured 
%                   along the blade, that is, the distance between the 
%                   nodes.
% foilwt          : Nfoil-by-Nel table.  Weights to use when computing
%                   airfoil coefficients from splined tables.
% aoaz            : Zero-lift angles-of-attack for each element. Should
%                   be computed precisely from the airfoil tables.
% aoast           : 2*Nel vector, containing deep-stall angles-of-attack
%                   for each element.  Alternating positive, negative.
% xas,yas         : X^a and Y^a coordinates of the reference section
%                   coordinate system.
% iq              : DOF indices for qy, qB, qn1, qn2 associated with each
%                   element.

[idofs,idofm,inods,inodm,Ndof] = getDOFRefs (s);

Nel = 0;
for ib = 1:s.Nb
   Nel = Nel + s.blade(ib).Nel;
end

Tas     = zeros(3,3*Nel);
ch      = zeros(Nel,1);
Lel     = zeros(Nel,1);
foilwt  = zeros(size(a.foilwt,1),Nel);
aoaz    = zeros(Nel,1);
aoast   = zeros(2*Nel,1);
xas     = zeros(Nel,1);
yas     = zeros(Nel,1);
iq      = zeros(24,Nel);
for ibod = 5:7

   jb3 = 3*(ibod-1);

   if (ibod == 5)
      idref = idofs(6);
      Neb   = s.blade(1).Nel;
      conns = s.blade(1).conn;
      elref = 0;
   elseif (ibod == 6)
      idref = idofs(7);
      Neb   = s.blade(2).Nel;
      conns = s.blade(2).conn;
      elref = s.blade(1).Nel;
   elseif (ibod == 7)
      idref = idofs(8);
      Neb   = s.blade(3).Nel;
      conns = s.blade(3).conn;
      elref = s.blade(1).Nel + s.blade(2).Nel;
   end

   for iel = 1:Neb

      jel = elref + iel;
      jc3   =  3*(jel-1);
      jc2   =  2*(jel-1);

      conn  = conns(:,iel);
      rdof  = idref + 6*(conn(1)-1);
      n1dof = idref + 6*(conn(2)-1);
      n2dof = idref + 6*(conn(3)-1);

      iq(1:6,jel)   = idofs(3)+[1:6];
      iq(7:12,jel)  = idref+[1:6];
      iq(13:18,jel) = n1dof+[1:6];
      iq(19:24,jel) = n2dof+[1:6];

      Tas(:,jc3+[1:3]) = a.Ta_s(:,3*(jel-1)+[1:3]);

      ch(jel) = a.chord(jel);
      Lel(jel) = a.Lel(jel);
      foilwt(:,jel) = a.foilwt(:,jel);
      aoast(jc2+[1:2]) = [0.5;-0.5];  % (Deactivated deep-stall tau anyways...)
      aoaz(jel) = aoazero (a.aoas,a.kfoils,foilwt(:,jel));

      xas(jel)  = (a.xpc(jel) - a.acent(jel))*ch(jel);
      yas(jel)  = 0;   % (For now...)

   end

end
