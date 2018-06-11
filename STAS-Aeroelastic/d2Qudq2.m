function d2Q = d2Qudq2 (qn,qB,Pn,PB)
%
% Compute the second partial derivative of the Qu matrix with respect to
% nodal and body reference degrees-of-freedom, for a single node.
%
% [v;w] = Qu dq/dt.
%
% Version:        Changes:
% --------        -------------
% 05.11.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 05.11.2017      Verified with complex step using dQudq.m.
% ERRORS FOUND FOR SOME SIMPLE INPUTS
%
% Inputs:
% -------
% qn              : DOFs associated with the node.
% qB              : DOFs associated with the body reference node.
% Pn              : Undeformed pos,rot of the node wrt body reference.
% PB              : Undeformed pos,rot of the body reference.
%
% Outputs:
% --------
% d2Q             : Qu is of dimension 6-by-12, and there are 12-by-12
%                   second derivatives.
%

d2Q = zeros(6,12*12*12);

r = qn(1:3) + Pn(1:3);

[Fr,TBB0,dTBB0,TB0g] = FRefMatrix (qB(4:6),PB(4:6));
d2TBB0 = d2Tdth2 (qB(4:6),TBB0,dTBB0);
d3TBB0 = d3Tdth3 (qB(4:6),TBB0,dTBB0,d2TBB0);
Fmat = Dmatrix (Fr);

TBg = TB0g*TBB0;
TgB = TBg.';
Tn0B = TFromTheta (Pn(4:6));

[F,Tnn0,dTnn0] = Fmatrix (qn(4:6),TBg,Tn0B);
d2Tnn0 = d2Tdth2 (qn(4:6),Tnn0,dTnn0);
Gmat = Dmatrix (F);

dFr = dFrefdth (qB(4:6),TBB0,dTBB0,d2TBB0,TB0g);
dFmat = zeros(3,3*3);
for jj = 1:3
   jc3 = 3*(jj-1);
   jc9 = 9*(jj-1);
   dFmat(:,jc3+[1:3]) = Dmatrix (dFr(:,jc9+[1:9]));
end

dF = dFdth (Tnn0,dTnn0,d2Tnn0,TBB0,dTBB0,TB0g,Tn0B);
dGmat = zeros(3,6*3);
for jj = 1:6
   jc3 = 3*(jj-1);
   jc9 = 9*(jj-1);
   dGmat(:,jc3+[1:3]) = Dmatrix (dF(:,jc9+[1:9]));   
end

d2Fr = d2Frefdth2 (qB(4:6),TBB0,dTBB0,d2TBB0,d3TBB0,TB0g);
d2Fmat = zeros(3,3*3*3);
for jj = 1:9
   jc3 = 3*(jj-1);
   jc9 = 9*(jj-1);
   d2Fmat(:,jc3+[1:3]) = Dmatrix (d2Fr(:,jc9+[1:9]));
end

d2F = d2Fdth2 (qB(4:6),qn(4:6),TB0g,Tn0B);
d2Gmat = zeros(3,3*6*6);
for jj = 1:36
   jc3 = 3*(jj-1);
   jc9 = 9*(jj-1);
   d2Gmat(:,jc3+[1:3]) = Dmatrix (d2F(:,jc9+[1:9]));
end

% d2Qi/ddj dPhik.  
% Same as d2Qi/dOj dPhik, but different jd index.
jd = 7;   % d
kd = 4;   % Phi

id = 4;   % Phi
ic = 12*12*(jd-1) + 12*(kd-1) + id-1;

for jj = 1:3

   jc144 = 12*12*(jj-1);

   for kk = 1:3

      kc12 = 12*(kk-1);
      kc3  =  3*(kk-1);

      for ii = 1:3

         ic3 = 3*(ii-1);
         ic9 = 9*(ii-1);

         d2Q(1:3,ic+jc144+kc12+ii) = (dTBB0(:,kc3+[1:3]).')*dTBB0(:,ic3+jj) ...
                                   + (TBB0.')*d2TBB0(:,ic9+kc3+jj);

      end

   end

end

% d2Qi/dPhij ddk.  
jd = 4;   % Phi
kd = 7;   % d

id = 4;   % Phi
ic = 12*12*(jd-1) + 12*(kd-1) + id-1;

for jj = 1:3

   jc144 = 12*12*(jj-1);
   jc3   =     3*(jj-1);

   for kk = 1:3

      kc12 = 12*(kk-1);

      for ii = 1:3

         ic3 = 3*(ii-1);
         ic9 = 9*(ii-1);

         d2Q(1:3,ic+jc144+kc12+ii) = (dTBB0(:,jc3+[1:3]).')*dTBB0(:,ic3+kk) ...
                                   + (TBB0.')*d2TBB0(:,ic9+jc3+kk);

      end

   end

end

% d2Qi/dPhij dPhik
jd = 4; % Phi
kd = 4; % Phi

for jj = 1:3

   jc144 = 144*(jj-1);
   jc27  =  27*(jj-1);
   jc18  =  18*(jj-1);
   jc9   =   9*(jj-1);
   jc6   =   6*(jj-1);
   jc3   =   3*(jj-1);

   dTj = dTBB0(:,jc3+[1:3]);
   dFj = dFmat(:,jc3+[1:3]);
   dGj = dGmat(:,jc3+[1:3]);

   for kk = 1:3

      kc12 = 12*(kk-1);
      kc9  =  9*(kk-1);
      kc6  =  6*(kk-1);
      kc3  =  3*(kk-1);

      dTk   = dTBB0(:,kc3+[1:3]);
      d2Tjk = d2TBB0(:,jc9+kc3+[1:3]);

      dFk   = dFmat(:,kc3+[1:3]);
      d2Fjk = d2Fmat(:,jc9+kc3+[1:3]);

      dGk   = dGmat(:,kc3+[1:3]);
      d2Gjk = d2Gmat(:,jc18+kc3+[1:3]);

      id = 1;   % Start with the columns associated with O.
      ic = 12*12*(jd-1) + 12*(kd-1) + id-1;

      % Rows associated with v.
      % (Ahat*That)_i = (A*That)_i = I 
      term1 = ((TB0g*d2Tjk).');  % TB0g: the tilde transforms are I, so the
      term2 = zeros(3);          % outer TB0g is not cancelled.
      term3 = zeros(3);
      term4 = zeros(3);
      term5 = zeros(3);
      term6 = zeros(3);

      d2Q(1:3,ic+jc144+kc12+[1:3]) = term1;

      % Rows associated with w are zero, as w is independent of O.

      id = 4;   % Columns associated with Phi.
      ic = 12*12*(jd-1) + 12*(kd-1) + id-1;

      % Rows associated with v.
      % (Ahat*That)_i = 0.  (A*That)_i = 0.
      term1 = zeros(3);
      term2 = zeros(3);
%     term3 = zeros(3);
%     term4 = zeros(3);
      term5 = zeros(3);
      term6 = zeros(3);
      for ii = 1:3
         ic3         = 3*(ii-1);
         dTi         = dTBB0(:,ic3+[1:3]);
         d2Tik       = d2TBB0(:,kc9+ic3+[1:3]);
         d2Tij       = d2TBB0(:,jc9+ic3+[1:3]);
         d3Tijk      = d3TBB0(:,jc27+kc9+ic3+[1:3]);
         term1(:,ii) = (d2Tjk.')* dTi  *r;
         term2(:,ii) =   (dTj.')*d2Tik *r;
         term5(:,ii) =   (dTk.')*d2Tij *r;
         term6(:,ii) =  (TBB0.')*d3Tijk*r;
      end

      d2Q(1:3,ic+jc144+kc12+[1:3]) = term1 + term2 + term5 + term6;

      % Rows associated with w.
      % (Ahat*That)_i = F*I.  A*That = 0.  
      term1 = ((TB0g*d2Tjk).')*Fmat;
      term2 = ((TB0g*dTj).')*dFk;
      term3 = ((TB0g*dTk).')*dFj;
      term4 = TgB*d2Fjk;
%     term5 = zeros(3);
%     term6 = zeros(3);

      d2Q(4:6,ic+jc144+kc12+[1:3]) = term1 + term2 + term3 + term4;

      id = 7;   % Columns associated with d.
      ic = 12*12*(jd-1) + 12*(kd-1) + id-1;

      % Rows associated with v.
      % (Ahat*That)_i = (A*That)_i = TBg.
      term1 = (d2Tjk.')*TBB0;
      term2 = (dTj.')*dTk;
      term3 = (dTk.')*dTj;
      term4 = (TBB0.')*d2Tjk;
%     term5 = zeros(3);
%     term6 = zeros(3);

      d2Q(1:3,ic+jc144+kc12+[1:3]) = term1 + term2 + term3 + term4;

      % Rows associated with w.  Zero, both Ahat and A contribute zero here.

      id = 10;  % Columns associated with theta.
      ic = 12*12*(jd-1) + 12*(kd-1) + id-1;

      % Rows associated with v.  
      % (Ahat*That)_i = 0.  (A*That)_i = 0.  Then, the terms multiplying
      % r are also zero; these depend on Phi, not theta.

      % Rows associated with w.
      % (Ahat*That)_i = G*I.  A*That = 0.
      term1 = ((TB0g*d2Tjk).')*Gmat;
      term2 = ((TB0g*dTj).')*dGk;
      term3 = ((TB0g*dTk).')*dGj;
      term4 = TgB*d2Gjk;
%     term5 = zeros(3);
%     term6 = zeros(3);

      d2Q(1:3,ic+jc144+kc12+[1:3]) = term1 + term2 + term3 + term4;

   end

end

% d2Qi/dPhij dthk
jd = 4;  % Phi
kd = 10; % theta

for jj = 1:3

   jc144 = 144*(jj-1);
   jc18  =  18*(jj-1);
   jc3   =   3*(jj-1);

   dTj = dTBB0(:,jc3+[1:3]);
   dGj = dGmat(:,jc3+[1:3]);

   for kk = 1:3

      kc12 = 12*(kk-1);
      kc3a =  3*(kk+3-1);

      dGk = dGmat(:,kc3a+[1:3]);
      d2Gjk = d2Gmat(:,jc18+kc3a+[1:3]);

      %id = 1;   % Columns associated with O.
      %ic = 12*12*(jd-1) + 12*(kd-1) + id-1;

      % Rows associated with v.
      % (Ahat*That)_i = (A*That)_i = I.  Derivatives of Ahat are zero.

      % Rows associated with w are zero.

      %id = 4;   % Columns associated with Phi.
      %ic = 12*12*(jd-1) + 12*(kd-1) + id-1;

      % Rows associated with v.
      % (Ahat*That)_i = 0.  (A*That)_i = 0.

      % Rows associated with w.
      % (Ahat*That)_i = F*I.  A*That = 0.  
      % Derivatives of F wrt theta are zero.

      %id = 7;   % Columns associated with d.
      %ic = 12*12*(jd-1) + 12*(kd-1) + id-1;

      % Rows associated with v.
      % (Ahat*That)_i = (A*That)_i = TBg.  Derivatives of Ahat are zero.

      % Rows associated with w.  Zero, both Ahat and A contribute zero.

      id = 10;   % Columns associated with theta.
      ic = 12*12*(jd-1) + 12*(kd-1) + id-1;

      % Rows associated with v.  
      % (Ahat*That)_i = 0.  (A*That)_i = 0.

      % Rows associated with w.
      % (Ahat*That)_i = G*I.  A*That = 0.
      term1 = ((TB0g*dTj).')*dGk;
      term2 = TgB*d2Gjk;

      d2Q(4:6,ic+jc144+kc12+[1:3]) = term1 + term2;

   end

end

% d2Qi/dthj dPhik
jd = 10;  % theta
kd = 4;   % Phi

for jj = 1:3

   jc144 = 144*(jj-1);
   jc18a = 18*(jj+3-1);
   jc3a  =  3*(jj+3-1);

   dGj = dGmat(:,jc3a+[1:3]);

   for kk = 1:3

      kc12 = 12*(kk-1);
      kc3 = 3*(kk-1);

      dTk = dTBB0(:,kc3+[1:3]);
      d2Gjk = d2Gmat(:,jc18a+kc3+[1:3]);

      id = 10;   % Columns associated with theta.
      ic = 12*12*(jd-1) + 12*(kd-1) + id-1;

      % Rows associated with w.
      % (Ahat*That)_i = G*I.  A*That = 0.
      term1 = ((TB0g*dTk).')*dGj;
      term2 = TgB*d2Gjk;

      d2Q(4:6,ic+jc144+kc12+[1:3]) = term1 + term2;

   end

end

% d2Qi/dthj dthk
jd = 10;  % theta
kd = 10;  % theta

for jj = 1:3

   jc144 = 144*(jj-1);
   jc18a = 18*(jj+3-1);

   for kk = 1:3

      kc12 = 12*(kk-1);
      kc3a = 3*(kk+3-1);

      d2Gjk = d2Gmat(:,jc18a+kc3a+[1:3]);

      % Only columns associated with theta and rows associated with w
      % have a nonzero value, since dAh/dth is nonzero only in this
      % case.

      id = 10;  % Columns associated with theta.
      ic = 12*12*(jd-1) + 12*(kd-1) + id-1;

      % Rows associated with w.
      % (Ahat*That)_i = G*I.  A*That = 0.
      d2Q(4:6,ic+jc144+kc12+[1:3]) = TgB*d2Gjk;

   end

end

d2Q = sparse(d2Q);
