function [dxdt,A] = TCfun (x,ret,x0,u,c);

x1 = x0;
x1(ret) = x;

[dxfdt,yout,AA,BB,CC,DD,blydof,bludof] =                           ...
      turbineControl (1,x1,u,c.cpar,c.cpct,                        ...
                      c.KeTab,c.WVTab,c.WPTab,c.bminTab,c.KTables, ...
                      c.KFTab,c.KSTab,c.KSqTab,c.KpiTab,c.KiiTab);


dxdt = dxfdt(ret);
A = AA(ret,ret);