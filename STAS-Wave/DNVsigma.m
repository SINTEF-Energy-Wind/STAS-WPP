function sigma = DNVsigma (f,Tp)

sigma = zeros(size(f));
sigma(f <= 1./Tp) = 0.07;
sigma(f > 1./Tp) = 0.09;