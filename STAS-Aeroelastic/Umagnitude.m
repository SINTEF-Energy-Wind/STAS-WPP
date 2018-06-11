function [Umag,Dy] = Umagnitude (Ua)

Umag = sqrt(Ua(1)^2 + Ua(2)^2);
Dy = [Ua(1)/Umag Ua(2)/Umag];