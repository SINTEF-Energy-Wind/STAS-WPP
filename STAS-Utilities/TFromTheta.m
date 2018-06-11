function T = TFromTheta (th)

TTH = vecToSpin (th);
T = expm (TTH);
