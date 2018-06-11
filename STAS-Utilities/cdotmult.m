function [zR,zI] = cdotmult (xR,xI,yR,yI)

zR = xR.*yR - xI.*yI;
zI = xR.*yI + xI.*yR;
