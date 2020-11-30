%{
if (Vinf(1) < 10)
   eval(["save('-binary','Lpsi_V0" int2str(round(10*Vinf(1))) ".bin','Lpsi');"]);
   eval(["save('-binary','Apsi_V0" int2str(round(10*Vinf(1))) ".bin','Apsi');"]);
   eval(["save('-binary','Bpsi_V0" int2str(round(10*Vinf(1))) ".bin','Bpsi');"]);
   eval(["save('-binary','Cpsi_V0" int2str(round(10*Vinf(1))) ".bin','Cpsi');"]);
   eval(["save('-binary','Dpsi_V0" int2str(round(10*Vinf(1))) ".bin','Dpsi');"]);
   eval(["save('-binary','dret_V0" int2str(round(10*Vinf(1))) ".bin','dret');"]);
   eval(["save('-binary','shape_V0" int2str(round(10*Vinf(1))) ".bin','shape0');"]);
else
   eval(["save('-binary','Lpsi_V" int2str(round(10*Vinf(1))) ".bin','Lpsi');"]);
   eval(["save('-binary','Apsi_V" int2str(round(10*Vinf(1))) ".bin','Apsi');"]);
   eval(["save('-binary','Bpsi_V" int2str(round(10*Vinf(1))) ".bin','Bpsi');"]);
   eval(["save('-binary','Cpsi_V" int2str(round(10*Vinf(1))) ".bin','Cpsi');"]);
   eval(["save('-binary','Dpsi_V" int2str(round(10*Vinf(1))) ".bin','Dpsi');"]);
   eval(["save('-binary','dret_V" int2str(round(10*Vinf(1))) ".bin','dret');"]);
   eval(["save('-binary','shape_V" int2str(round(10*Vinf(1))) ".bin','shape0');"]);
end
%}


txt = outnm;

if (Vmag < 10)
   eval(["save('-ascii','xpsi" txt "V0" int2str(round(10*Vmag)) ".txt','xpsi');"]);
   eval(["save('-ascii','upsi" txt "V0" int2str(round(10*Vmag)) ".txt','upsi');"]);
   eval(["save('-ascii','dxpsi" txt "V0" int2str(round(10*Vmag)) ".txt','dxpsi');"]);
   eval(["save('-ascii','ypsi" txt "V0" int2str(round(10*Vmag)) ".txt','ypsi');"]);
   eval(["save('-ascii','Rgrav" txt "V0" int2str(round(10*Vmag)) ".txt','Rgrav');"]);
   eval(["save('-binary','Lpsi" txt "V0" int2str(round(10*Vmag)) ".bin','Lpsi');"]);
   eval(["save('-binary','Apsi" txt "V0" int2str(round(10*Vmag)) ".bin','Apsi');"]);
   eval(["save('-binary','Bpsi" txt "V0" int2str(round(10*Vmag)) ".bin','Bpsi');"]);
   eval(["save('-binary','Cpsi" txt "V0" int2str(round(10*Vmag)) ".bin','Cpsi');"]);
   eval(["save('-binary','Dpsi" txt "V0" int2str(round(10*Vmag)) ".bin','Dpsi');"]);
   eval(["save('-binary','dret" txt "V0" int2str(round(10*Vmag)) ".bin','dret');"]);
   eval(["save('-binary','shape" txt "V0" int2str(round(10*Vmag)) ".bin','shape0');"]);
   eval(["save('-binary','mdamp" txt "V0" int2str(round(10*Vmag)) ".bin','mdamp0');"]);
   eval(["save('-binary','bldof" txt "V0" int2str(round(10*Vmag)) ".bin','bldof');"]);
else
   eval(["save('-ascii','xpsi" txt "V" int2str(round(10*Vmag)) ".txt','xpsi');"]);
   eval(["save('-ascii','upsi" txt "V" int2str(round(10*Vmag)) ".txt','upsi');"]);
   eval(["save('-ascii','dxpsi" txt "V" int2str(round(10*Vmag)) ".txt','dxpsi');"]);
   eval(["save('-ascii','ypsi" txt "V" int2str(round(10*Vmag)) ".txt','ypsi');"]);
   eval(["save('-ascii','Rgrav" txt "V" int2str(round(10*Vmag)) ".txt','Rgrav');"]);
   eval(["save('-binary','Lpsi" txt "V" int2str(round(10*Vmag)) ".bin','Lpsi');"]);
   eval(["save('-binary','Apsi" txt "V" int2str(round(10*Vmag)) ".bin','Apsi');"]);
   eval(["save('-binary','Bpsi" txt "V" int2str(round(10*Vmag)) ".bin','Bpsi');"]);
   eval(["save('-binary','Cpsi" txt "V" int2str(round(10*Vmag)) ".bin','Cpsi');"]);
   eval(["save('-binary','Dpsi" txt "V" int2str(round(10*Vmag)) ".bin','Dpsi');"]);
   eval(["save('-binary','dret" txt "V" int2str(round(10*Vmag)) ".bin','dret');"]);
   eval(["save('-binary','shape" txt "V" int2str(round(10*Vmag)) ".bin','shape0');"]);
   eval(["save('-binary','mdamp" txt "V" int2str(round(10*Vmag)) ".bin','mdamp0');"]);
   eval(["save('-binary','bldof" txt "V" int2str(round(10*Vmag)) ".bin','bldof');"]);
end
