function [mu,dmu,d2mu,TsB,dTsB,d2TsB] = ...
                  prepareElementInput (linFlag,qB,qn1,qn2,PB,Pn1,Pn2);

[xe,TsB] = elementCSFromNodes (qn1,qn2,Pn1,Pn2);
dTsB = derivElementCS (qn1,qn2,Pn1,Pn2,TsB);
mu = getMu (qn1,qn2,Pn1,Pn2,TsB);
dmu = dmudq (qn1,qn2,Pn1,Pn2,TsB,dTsB);

if (linFlag == 1)
   d2TsB = secondDerivElementCS (qn1,qn2,Pn1,Pn2,TsB);
   d2mu = d2mudq2 (qn1,qn2,Pn1,Pn2,TsB,dTsB,d2TsB);
else
   d2TsB = sparse(3,3*12*12);
   d2mu = sparse(12,12*12);
end
