maxTotEnergyVect = zeros(1000,1);
for i = 1:1000
    maxTotEnergyVect(i) = WaveEquation1D(2000+80*i, 1, 100);
end