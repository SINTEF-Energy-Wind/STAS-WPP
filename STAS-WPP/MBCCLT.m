function [Lpsi,Rpsi,ypsi,Apsi,Bpsi,Cpsi,Dpsi] =                ...
                  MBCCLT (xpsi,upsi,Neta,                      ...
                          TpsixB,dTpsixB,TBxpsi,TpsiuB,TBypsi, ...
                          s,a,epar,ppar,ypar,c,                ...
                          linFlag,psiFlag,modFlag,shpFlag,     ...                                      
                          grav,P,ret,slv,shape0,mdamp0,        ...
                          Tas,Try,ch,Lel,foilwt,aoaz,aoast,    ...
                          xas,yas,Psi0,igen,ipit,iyaw)
%
% Return the multi-blade coordinate transformed version of the dynamic
% closed-loop turbine equations.
%
% Version:        Changes:
% --------        -------------
% 14.02.2019      Original code.
%
% Version:        Verification:
% --------        -------------
% 14.02.2019      
%

% Transform the input xpsi and upsi to body coordinates.
x = TpsixB*xpsi;
u = TpsiuB*upsi;

% Call the CLT function in body coordinates.
[L,R,yout,A,B,C,D,jnkblx,jnkblu,jnkbly] =                             ...
      buildClosedLoopTurbine (x,u,s,a,epar,ppar,ypar,c,               ...
                              linFlag,psiFlag,modFlag,shpFlag,        ...                                      
                              grav,P,ret,slv,shape0,mdamp0,           ...
                              Tas,Try,ch,Lel,foilwt,aoaz,aoast,       ...
                              xas,yas,Psi0,igen,ipit,iyaw);

% Get the rotor speed for the MBC transform.  This is by definition 
% the time derivative of the DOF used to define the azimuth.  Here
% this is assumed to be the generator rotor.  [Could also be the hub
% center node.]
iW = 2*Neta - 4;
W  = x(iW);

Lpsi = TBxpsi*L*TpsixB;
Rpsi = TBxpsi*(R - W*L*dTpsixB*xpsi);
ypsi = TBypsi*yout;

if (linFlag == 1)

   LdT  = L*dTpsixB;
   Apsi = TBxpsi*(A*TpsixB - W*LdT);
   Apsi(:,iW) = Apsi(:,iW) - TBxpsi*LdT*xpsi;
   Bpsi = TBxpsi*B*TpsiuB;
   Cpsi = TBypsi*C*TpsixB;
   Dpsi = TBypsi*D*TpsiuB;

else

   Apsi = sparse(size(A,1),size(A,2));
   Bpsi = sparse(size(B,1),size(B,2));
   Cpsi = sparse(size(C,1),size(C,2));
   Dpsi = sparse(size(D,1),size(D,2));

end
