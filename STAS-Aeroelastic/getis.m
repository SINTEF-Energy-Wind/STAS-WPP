function [iad,ia1,ia2,iVih,iVi,               ...
    iq,idq,ixng1,ixng2,ivng1,ivng2,           ...
    iwg,iwa,iaq,iC,iFldm,iFa,iFp,iFr,iFzts,   ...
    iVg,iUa,iUmag,iUr,iUzts,iVzts,iViq,       ...
    iViy,iVixyz,iWmag,ixeg,ixhg,              ...
    ixnr1,ixnr2,ixer,ir,iL,iz,iPr] = getis (iel);
% Version:        Changes:
% --------        -------------
% 10.02.2018      Original code.
%
% Version:        Verification:
% --------        -------------
% 10.02.2018      Checked visually.

ioff = 7*(iel-1);

iad   = ioff;
ia1   = iad + 1;
ia2   = ia1 + 1;
iVih  = ia2 + 1;
iVi   = iVih + 2;

ioff = 137*(iel-1);

iq    = ioff;
idq   = iq    + 24;
ixng1 = idq   + 24; 
ixng2 = ixng1 + 3;
ivng1 = ixng2 + 3;
ivng2 = ivng1 + 3;
iwg   = ivng2 + 3;
iwa   = iwg   + 3;
iaq   = iwa   + 3;
iC    = iaq   + 1;
iFldm = iC    + 3; 
iFa   = iFldm + 3;
iFp   = iFa   + 6;
iFr   = iFp   + 6;
iFzts = iFr   + 6; 
iVg   = iFzts + 3;
iUa   = iVg   + 3;
iUmag = iUa   + 3; 
iUr   = iUmag + 1;
iUzts = iUr   + 3;
iVzts = iUzts + 3;
iViq  = iVzts + 3;
iViy  = iViq  + 2;
iVixyz= iViy  + 2;
iWmag = iVixyz+ 3;
ixeg  = iWmag + 1;
ixhg  = ixeg  + 3;
ixnr1 = ixhg  + 3;
ixnr2 = ixnr1 + 3;
ixer  = ixnr2 + 3;
ir    = ixer  + 3;
iL    = ir    + 1;
iz    = iL    + 1;
iPr   = iz    + 1;
