function [mu,dmu,d2mu,TsB,dTsB,d2TsB] = ...
                  prepareElementInput (qB,qn1,qn2,PB,Pn1,Pn2);

[xe,TsB] = elementCSFromNodes (qn1,qn2,Pn1,Pn2);
dTsB = derivElementCS (qn1,qn2,Pn1,Pn2,TsB);
d2TsB = secondDerivElementCS (qn1,qn2,Pn1,Pn2,TsB);

mu = getMu (qn1,qn2,Pn1,Pn2,TsB);
dmu = dmudq (qn1,qn2,Pn1,Pn2,TsB,dTsB);
d2mu = d2mudq2 (qn1,qn2,Pn1,Pn2,TsB,dTsB,d2TsB);
