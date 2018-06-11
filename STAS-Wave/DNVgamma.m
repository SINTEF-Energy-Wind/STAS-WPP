function gamma = DNVgamma (Hs,Tp)

TpsqHs = Tp./sqrt(Hs);
gamma = zeros(size(TpsqHs));
gamma(TpsqHs <= 3.6) = 5;
gamma(TpsqHs >= 5) = 1;
gamma((TpsqHs > 3.6) && (TpsqHs < 5)) = exp(5.75 - 1.15*TpsqHs);


