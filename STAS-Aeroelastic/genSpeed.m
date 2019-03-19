function [Wg,Dy] = genSpeed (qgen,dqgen,Pgen)
%
% Find the present generator speed Wg and the linearizations dWg/dq
% and dWg/dqdot.
%
%   States:           y vector:         u vector:
%                     qBs       1:6
%                     qns       7:12
%                     qBr      13:18
%                     qnr      19:24
%                     dq/dt    25:48
%
% Version:        Changes:
% --------        -------------
% 17.01.2019      Original code.
%
% Version:        Verification:
% --------        -------------
% 17.01.2019      Derivatives verified by complex step.
%
% Inputs:
% -------
% qgen              : [qBs,qns,qBr,qnr], the nodal and body reference
%                     DOFs for the generator stator and rotor. 
% dqgen             : Ditto for velocities.
% Pgen              : Ditto, nodal positions.
%
% Outputs:
% --------
% Wg                : Gen. shaft speed.
% Dy                : dWg/dqBs, dWg/dqns, dWg/dqBr, dWg/dqnr,
%                     dWg/dqBsd, dWg/dqnsd, dWg/dqBrd, dWg/dqnrd.

Dy = zeros(1,48);

[vgs,Dys] = globalVelocity (qgen(7:12),qgen(1:6),Pgen(7:12),Pgen(1:6), ...
                            dqgen(7:12),dqgen(1:6));
[vgr,Dyr] = globalVelocity (qgen(19:24),qgen(13:18),Pgen(19:24),Pgen(13:18), ...
                            dqgen(19:24),dqgen(13:18));

Td0g = TFromTheta (Pgen(16:18));
[Tdd0,dTdd0] = dTdth (qgen(16:18));
Tdg = Td0g*Tdd0;

Tvec = Tdg(:,3).';
vdiff = vgr(4:6) - vgs(4:6);

Wg = Tvec*vdiff;

Dy(1:6)   = -Tvec*Dys(4:6,1:6);
Dy(7:12)  = -Tvec*Dys(4:6,7:12);
Dy(25:30) = -Tvec*Dys(4:6,13:18);
Dy(31:36) = -Tvec*Dys(4:6,19:24);
Dy(13:18) = Tvec*Dyr(4:6,1:6);
Dy(19:24) = Tvec*Dyr(4:6,7:12);
Dy(37:42) = Tvec*Dyr(4:6,13:18);
Dy(43:48) = Tvec*Dyr(4:6,19:24);

for ith = 1:3
   ic3 = 3*(ith-1);
   Dy(3+ith) = Dy(3+ith) + ((Td0g*dTdd0(:,ic3+3)).')*vdiff;
end 
