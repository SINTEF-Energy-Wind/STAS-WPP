function iels = getiels (s)
% Version:        Changes:
% --------        -------------
% 01.11.2017      Original code.
%
% Version:        Verification:
% --------        -------------
% 01.11.2017      Checked visually.

iels(1) = 0;
iels(2) = s.Nef;
iels(3) = iels(2) + s.Net;
iels(4) = iels(3) + s.Nev + s.Ner + s.Nen;
iels(5) = 0;                            % The front bearing.  Keeps the
iels(6) = iels(4) + s.Ned + s.Neh + 3;  % indices consistent with inods.
iels(7) = iels(6) + s.Neb;
iels(8) = iels(7) + s.Neb;
