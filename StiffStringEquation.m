clear all;
close all;
clc;

fs = 44100;     % Sampling rate
E = 2E11;       % Young's Modulus
L = 1.5;       % String length
rho = 7850;     % Density of steel [kg/m^3]
gamma = 500;        % Speed
k = 1/fs;
s0 = 0.0;
s1 = 0.00;
numKappa = 10000;
kappaLoop = false;
if kappaLoop
    muSq = zeros(numKappa,1);
    lambdaSq = zeros(numKappa,1);
    for kappa = 1 : numKappa
        % kappa = 0.5; %sqrt((E * inertia) / rho * A * L^4);  % Stiffness
        % h = nthroot(k^2 * kappa^2 + 16 * gamma^2 * k^2, 4);
        h = sqrt((gamma^2*k^2 + 4 * s1 * k + sqrt((gamma^2 * k^2 + 4 * s1 * k)^2 + 16 * (kappa / 100)^2 * k^2)) / 2);
        muSq(kappa) = (k * (kappa / 100) / h^2)^2;
        lambdaSq(kappa) = (gamma * k / h)^2;
    end
    plot(lambdaSq)
    hold on; plot(muSq)
    plot(lambdaSq + 4*muSq)
else
     kappa = 2;
     h = sqrt((gamma^2*k^2 + 4 * s1 * k + sqrt((gamma^2 * k^2 + 4 * s1 * k)^2 + 16 * kappa^2 * k^2)) / 2);
     muSq = (k * kappa / h^2)^2;
     lambdaSq = (gamma * k / h)^2;
end
P = 1/2;

N = floor(L / h) + 1;

cosWidth = round(N / 15);
raisedCos = 0.5 * (cos(pi:(2*pi)/cosWidth:3*pi) + 1);
PIdx = floor (P * N);
u = zeros(N, 1);
u(ceil(N * P - cosWidth / 2 : ceil(N * P + cosWidth / 2))) = raisedCos;
% u(3:N-2) = rand(N - 4,1) - 0.5; 
uPrev = u;
plot (u);
uNext = zeros(N,1);

lengthSound = fs;
out = zeros(lengthSound, 1);

kinEnergy = zeros(lengthSound,1);
potEnergy = zeros(lengthSound,1);

drawString = false;

bc = "ss";

for t = 1 : lengthSound
    if bc == "clamped"
        for l = 2 : N - 1
            if l > 2 && l < N - 1
                  uNext(l) = (2 * u(l) - uPrev(l) + lambdaSq * (u(l+1) -2 * u(l) + u(l-1)) ...
                      - muSq * (u(l+2) - 4*u(l+1) + 6*u(l) - 4*u(l-1) + u(l-2)) ...
                      + s0 * k * uPrev(l) + (2*s1*k)/h^2 ...
                      * ((u(l+1) - 2 * u(l) + u(l-1)) ...
                      - (uPrev(l+1) - 2 * uPrev(l) + uPrev(l-1)))) ...
                      / (1 + s0 * k);
            end
            potEnergy(t) = potEnergy(t) ...
                + gamma^2 / 2 * sum (1 / h * (u(l + 1) - u(l)) * (uPrev(l + 1) - uPrev(l))) ...
                + kappa^2 / 2 * 1/h^3 * sum((u(l + 1) - 2 * u(l) + u(l - 1)) ...
                * (uPrev(l + 1) - 2 * uPrev(l) + uPrev(l - 1)));
            kinEnergy(t) = kinEnergy(t) + 1 / 2 * sum (h * ((1 / k * (u(l) - uPrev(l)))^2));

        end
    elseif bc == "ss"
        for l = 1 : N
            if l > 2 && l < N - 1
                uNext(l) = (2 * u(l) - uPrev(l) + lambdaSq * (u(l+1) -2 * u(l) + u(l-1)) ...
                      - muSq * (u(l+2) - 4*u(l+1) + 6*u(l) - 4*u(l-1) + u(l-2)) ...
                      + s0 * k * uPrev(l) + (2*s1*k)/h^2 ...
                      * ((u(l+1) - 2 * u(l) + u(l-1)) ...
                      - (uPrev(l+1) - 2 * uPrev(l) + uPrev(l-1)))) ...
                      / (1 + s0 * k);
            elseif l == 2
                  uNext(l) = (2 * u(l) - uPrev(l) + lambdaSq * (u(l+1) -2 * u(l) + u(l-1)) ...
                      - muSq * (u(l+2) - 4*u(l+1) + 5*u(l) - 4*u(l-1)) ...
                      + s0 * k * uPrev(l) + (2*s1*k)/h^2 ...
                      * ((u(l+1) - 2 * u(l) + u(l-1)) ...
                      - (uPrev(l+1) - 2 * uPrev(l) + uPrev(l-1)))) ...
                      / (1 + s0 * k);
            elseif l == N - 1
                   uNext(l) = (2 * u(l) - uPrev(l) + lambdaSq * (u(l+1) -2 * u(l) + u(l-1)) ...
                      - muSq * (-4*u(l+1) + 5*u(l) - 4*u(l-1) + u(l-2)) ...
                      + s0 * k * uPrev(l) + (2*s1*k)/h^2 ...
                      * ((u(l+1) - 2 * u(l) + u(l-1)) ...
                      - (uPrev(l+1) - 2 * uPrev(l) + uPrev(l-1)))) ...
                      / (1 + s0 * k);
            end
            if l > 1 && l < N
                waveEqPotEnergy = gamma^2 / 2 * 1 / h * ...
                    (u(l) - u(l-1)) .* (uPrev(l) - uPrev(l-1));
                barPotEnergy = kappa^2 / 2 * 1/h^3 ...
                        * (u(l+1) - 2 * u(l) + u(l-1)) * (uPrev(l+1) - 2 * uPrev(l) + uPrev(l-1));
                potEnergy(t) = potEnergy(t) + waveEqPotEnergy + barPotEnergy;
            end
            if l == 1
                waveEqPotEnergy = 0;
                barPotEnergy = kappa^2 / 2 * 1/h^3 ...
                        * (-2 * u(1)) * (-2 * uPrev(1));
                potEnergy(t) = potEnergy(t) + waveEqPotEnergy + barPotEnergy;
            elseif l == N
                waveEqPotEnergy = gamma^2 / 2 * 1 / h * ...
                    (u(l) - u(l-1)) .* (uPrev(l) - uPrev(l-1));
                barPotEnergy = kappa^2 / 2 * 1/h^3 ...
                        * (-2 * u(N)) * (-2 * uPrev(N));
                potEnergy(t) = potEnergy(t) + waveEqPotEnergy + barPotEnergy;
            end
            kinEnergy(t) = kinEnergy(t) + 1 / 2 * h * ((1 / k * (u(l) - uPrev(l)))^2);
        end
    elseif bc == "free"
        for l = 1 : N
            if l > 2 && l < N - 1
                  uNext(l) = (2 * u(l) - uPrev(l) + lambdaSq * (u(l+1) -2 * u(l) + u(l-1)) ...
                      - muSq * (u(l+2) - 4*u(l+1) + 6*u(l) - 4*u(l-1) + u(l-2)) ...
                      + s0 * k * uPrev(l) + (2*s1*k)/h^2 ...
                      * ((u(l+1) - 2 * u(l) + u(l-1)) ...
                      - (uPrev(l+1) - 2 * uPrev(l) + uPrev(l-1)))) ...
                      / (1 + s0 * k);
            elseif l == 2
                  uNext(l) = (2 * u(l) - uPrev(l) + lambdaSq * (u(l+1) -2 * u(l) + u(l-1)) ...
                      - muSq * (u(l+2) - 4*u(l+1) + 7*u(l) - 4*u(l-1)) ...
                      + s0 * k * uPrev(l) + (2*s1*k)/h^2 ...
                      * ((u(l+1) - 2 * u(l) + u(l-1)) ...
                      - (uPrev(l+1) - 2 * uPrev(l) + uPrev(l-1)))) ...
                      / (1 + s0 * k);
            elseif l == N - 1
                   uNext(l) = (2 * u(l) - uPrev(l) + lambdaSq * (u(l+1) -2 * u(l) + u(l-1)) ...
                      - muSq * (-4*u(l+1) + 7*u(l) - 4*u(l-1) + u(l-2)) ...
                      + s0 * k * uPrev(l) + (2*s1*k)/h^2 ...
                      * ((u(l+1) - 2 * u(l) + u(l-1)) ...
                      - (uPrev(l+1) - 2 * uPrev(l) + uPrev(l-1)))) ...
                      / (1 + s0 * k);
            elseif l == 1
                   uNext(l) = (2 * u(l) - uPrev(l) + lambdaSq * (2 * u(l+1) - 2 * u(l)) ...
                      - muSq * (2 * u(l+2) - 8*u(l+1) + 6*u(l)) ...
                      + s0 * k * uPrev(l) + (2*s1*k)/h^2 ...
                      * ((2 * u(l+1) - 2 * u(l)) ...
                      - (2 * uPrev(l+1) - 2 * uPrev(l)))) ...
                      / (1 + s0 * k);
            elseif l == N
                   uNext(l) = (2 * u(l) - uPrev(l) + lambdaSq * (-2 * u(l) + 2 * u(l-1)) ...
                      - muSq * (6*u(l) - 8 * u(l-1) + 2 * u(l-2)) ...
                      + s0 * k * uPrev(l) + (2*s1*k)/h^2 ...
                      * ((-2 * u(l) + 2 * u(l-1)) ...
                      - (-2 * uPrev(l) + 2 * uPrev(l-1)))) ...
                      / (1 + s0 * k);
            end
            if l > 1 && l < N
                potEnergy(t) = potEnergy(t) ...
                    + gamma^2 / 2 * sum (1 / h * (u(l + 1) - u(l)) * (uPrev(l + 1) - uPrev(l))) ...
                    + kappa^2 / 2 * 1/h^3 * sum((u(l + 1) - 2 * u(l) + u(l - 1)) ...
                    * (uPrev(l + 1) - 2 * uPrev(l) + uPrev(l - 1)));
                kinEnergy(t) = kinEnergy(t) + 1 / 2 * sum (h * ((1 / k * (u(l) - uPrev(l)))^2));
            elseif l == 1
                potEnergy(t) = potEnergy(t) ...
                    + gamma^2 / 2 * sum (1 / h * (u(l + 1) - u(l)) * (uPrev(l + 1) - uPrev(l))) ...
                    + kappa^2 / 2 * 1/h^3 * sum((2 * u(l + 1) - 2 * u(l)) ...
                    * (2 * uPrev(l + 1) - 2 * uPrev(l)));
                kinEnergy(t) = kinEnergy(t) + 1 / 4 * sum (h * ((1 / k * (u(l) - uPrev(l)))^2));
            elseif l == N
                potEnergy(t) = potEnergy(t) ...
                    + gamma^2 / 2 * sum (1 / h * (u(l - 1) - u(l)) * (uPrev(l - 1) - uPrev(l))) ...
                    + kappa^2 / 2 * 1/h^3 * sum((- 2 * u(l) + 2 * u(l - 1)) ...
                    * (- 2 * uPrev(l) + 2 * uPrev(l - 1)));
                kinEnergy(t) = kinEnergy(t) + 1 / 4 * sum (h * ((1 / k * (u(l) - uPrev(l)))^2));
            end
        end
    end
%     plot3(cos([1:N]*2*pi/(N-1)), sin([1:N]*2*pi/(N-1)), u);
%     grid on
%     zlim([-1, 1]);
% totEnergy = kinEnergy + potEnergy;
% totEnergy = (totEnergy-totEnergy(1))/totEnergy(1);
% clf
% plot(kinEnergy(1:t)); hold on; plot(potEnergy(1:t)); plot(totEnergy(1:t)*10e5);
%     drawnow;
    out(t) = u(floor(N*P));
    uPrev = u;
    u = uNext;
end
totEnergy = kinEnergy + potEnergy;
totEnergy = (totEnergy-totEnergy(1))/totEnergy(1);
plot(totEnergy);
winSize = 2^12;
overlap = 2^8;
spectrogram(out*100,hanning(winSize),overlap,[],fs, 'MinThresh', -100,'yaxis')
% ylim([0,fs/2])
% plot(out);