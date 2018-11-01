clear all;
close all;
clc;

fs = 44100 * 2;     % Sampling rate
E = 2E11;       % Young's Modulus
L = 2;       % String length
rho = 7850;     % Density of steel [kg/m^3]
gamma = 300;        % Speed
k = 1/fs;
s0 = 0.0;
s1 = 0.0;

drawString = true;
flag = true;
kappa = 0;
h = sqrt((gamma^2*k^2 + 4 * s1 * k + sqrt((gamma^2 * k^2 + 4 * s1 * k)^2 + 16 * kappa^2 * k^2)) / 2);
gamma = gamma;
kappa = sqrt(((2*h^2-gamma^2*k^2-4*s1*k)^2-(gamma^2*k^2+4*s1*k)^2)/(16*k^2));
muSq = (k * kappa / h^2)^2;
lambdaSq = (gamma * k / h)^2;


P = 9/10;
N = floor(L / h) + 1;
cosWidth = round(N / 15);
raisedCos = 0.5 * (cos(pi:(2*pi)/cosWidth:3*pi) + 1);
PIdx = floor (P * N);
u = zeros(N, 1);
u(ceil(N * P - cosWidth / 2 : ceil(N * P + cosWidth / 2))) = raisedCos;
% u(3:N-2) = rand(N - 4,1) - 0.5; 
% u(N/2 : N/2+1) = 1;
uPrev = u;
uNext = zeros(N,1);

lengthSound = fs;
out = zeros(lengthSound, 1);

kinEnergy = zeros(lengthSound,1);
potEnergy = zeros(lengthSound,1);

drawString = true;
bc = "clamped";

mus=.4;                %static friction coeff
mud=.02;                %dynamic friction coeff (must be < mus!!)

nf = 0.001;

% sf=mus*nf;             % stiction force
% cf=mud*nf;             % coulomb force
vb = 0.2;
fb = 50;
% minVel = 0;
% maxDisp = 0.5;
% mode = "stick";
initBowPoint = P*N;
% bowArea = bowPoint:bowPoint+3;

%%
bowPoint = 30;
for t = 1 : lengthSound
    if bc == "clamped"
%         bowPoint = floor(5 + (N - 10) * (cos(4/(lengthSound/fs) * pi * t/fs + pi) + 1) * 0.5);
%         (u(bowPoint) - uPrev(bowPoint))
        fb = 0.03;%3 * t / lengthSound;

        vrel(t) =  1/k * (u(bowPoint) - uPrev(bowPoint)) - vb;
        phiRes(t) = phi(vrel(t));
        
%         if mode == "stick"
%             if u(bowPoint) >= maxDisp 
%                 mode = "slip";
%             else
%                 u(bowArea) = u(bowArea) + [vb/2 vb vb vb/2]';
%             end
%         end
%         if mode == "slip" && vel(t) >= minVel
%             mode = "stick";
%         end
        for l = 2 : N - 1
            if l > 2 && l < N - 1
                  uNext(l) = (2 * u(l) - uPrev(l) + lambdaSq * (u(l+1) -2 * u(l) + u(l-1)) ...
                      - muSq * (u(l+2) - 4*u(l+1) + 6*u(l) - 4*u(l-1) + u(l-2)) ...
                      + s0 * k * uPrev(l) + (2*s1*k)/h^2 ...
                      * ((u(l+1) - 2 * u(l) + u(l-1)) ...
                      - (uPrev(l+1) - 2 * uPrev(l) + uPrev(l-1)))) ...
                      / (1 + s0 * k);
%                   if l == bowPoint
%                        uNext(l) = uNext(l) + phiRes(t) * fb;
%                   end
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
            
            % Energy analysis
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
                waveEqPotEnergy = gamma^2 / 2 * 1 / h * ...
                    (u(l) - u(l-1)) .* (uPrev(l) - uPrev(l-1));
                barPotEnergy = kappa^2 / 2 * 1/h^3 ...
                        * (u(l+1) - 2 * u(l) + u(l-1)) * (uPrev(l+1) - 2 * uPrev(l) + uPrev(l-1));
                potEnergy(t) = potEnergy(t) + waveEqPotEnergy + barPotEnergy;
                kinEnergy(t) = kinEnergy(t) + 1 / 2 * h * ((1 / k * (u(l) - uPrev(l)))^2);
            end
            if l == 1
%                 waveEqPotEnergy = gamma^2 / 2 * 1 / h * ...
%                     (u(l) - u(l+1)) .* (uPrev(l) - uPrev(l+1));
                waveEqPotEnergy = 0;
                barPotEnergy = kappa^2 / 2 * 1/h^3 ...
                        * 0.5 * (2 * u(l+1) - 2 * u(l)) * (2 * uPrev(l+1) - 2 * uPrev(l));
                potEnergy(t) = potEnergy(t) + waveEqPotEnergy + barPotEnergy;
                kinEnergy(t) = kinEnergy(t) + 1 / 2 * 1 / 2 * h * ((1 / k * (u(l) - uPrev(l)))^2);
        
            elseif l == N
                waveEqPotEnergy = gamma^2 / 2 * 1 / h * ...
                    (u(l) - u(l-1)) .* (uPrev(l) - uPrev(l-1));
                barPotEnergy = kappa^2 / 2 * 1/h^3 ...
                        * 0.5 * (-2 * u(l) + 2 * u(l-1)) * (-2 * uPrev(l) + 2 * uPrev(l-1));
                potEnergy(t) = potEnergy(t) + waveEqPotEnergy + barPotEnergy;
                kinEnergy(t) = kinEnergy(t) + 1 / 2 * 1 / 2 * h * ((1 / k * (u(l) - uPrev(l)))^2);
        
            end
        end
    end
%     plot3(cos([1:N]*2*pi/(N-1)), sin([1:N]*2*pi/(N-1)), u);
%     grid on
%     zlim([-1, 1]);
% if flag == true
%     fret(t) = (sin(2*pi*0.5*t/fs + 0.5*pi) * 0.5 + 0.5) * -12 - 2;
%     if round(fret(t)) == -4
%         flag = false;
%     end
% else
%     fret(t) = -4;
% end

% if t > 0 && t <= 10000
%     fret = -12;
% elseif t > 10000 && t <= 20000
%     fret = -14;
% elseif t > 20000 
%     fret = -12;
% end
fret = -12;
fretPos = N - nthroot(2, 12)^fret * N;
% fretPos = 156;
pickup = floor(N/3);
% pickup = 110;
if drawString == true && mod(t,10) == 0 %&& t > 2*fs
    subplot(2,1,1)
    
    plot(u, 'Linewidth', 1);
%     ylim([-2e-19 2e-19])
    subplot(2,1,2)
    plot(phiRes(1:t))
%     pickupScat = [pickup, 0; 290, 0.67];
%     fretposScat = [fretPos, 0; 290, 0.60];
%     hold on; scatter(pickupScat(:,1), pickupScat(:,2));
%     scatter(fretposScat(:,1), fretposScat(:,2));
%     text (300, 0.67, "Pickup", 'Fontsize', 14);
%     text (300, 0.60, "Finger", 'Fontsize', 14);
% %     legend('Pick-up', 'Finger');
%     ylim([-0.75 0.75])
    drawnow;
end

fingerwidth = 4;
% if fret - fingerwidth > 0
fracPos = fretPos - floor(fretPos);
damping = 0;
% if fretPos ~= 0
%     for i = -fingerwidth/2:fingerwidth/2
%         uNext(floor(fretPos) + i) = damping * (uNext(floor(fretPos + i)) * (1-fracPos) + uNext(floor(fretPos) + 1 + i) * fracPos);
%     end
% end
%     uNext(floor(fretPos)) = damping * (uNext(floor(fretPos)) * (1-fracPos) + uNext(floor(fretPos) + 1) * fracPos);
% uNext(floor(fretPos)) = damping * uNext(floor(fretPos));
%     pos = sin(2*pi*t*2/fs) * 30 + 40;
%     fracPos = pos - floor(pos);
%     out(t) = uNext(floor(pos)) * (1-fracPos) + uNext(floor(pos)+1) * fracPos;
    out(t) = uNext(bowPoint);
    uPrev = u;
    u = uNext;
end
totEnergy = kinEnergy + potEnergy;
totEnergy = (totEnergy-totEnergy(1))/totEnergy(1);
winSize = 2^12;
overlap = 2^8;
spectrogram(out*100,hanning(winSize),overlap,[],fs, 'MinThresh', -100,'yaxis')
% ylim([0,fs/2])
plot(out);
% plot(vrel);

function [phiRes] = phi(eta)
    a = 100;
    epsilon = 0;
    muD = 0.3;
    muS = 0.8;
    v0 = 0.2;
%     phiRes = sign(eta) .* (epsilon + (1 - epsilon) * exp(-a*abs(eta)));
%     phiRes = sqrt(2*a) * eta * exp(-a*eta^2 + 1/2);

%         phiRes = 0.4 * exp(-abs(eta) / 0.01) + 0.45 * exp(-abs(eta) / 0.1) + 0.35; %stefania model 1
%     phiRes = muD + ((muS - muD) * v0) / (v0 + eta); %stefania model 2
    phiRes = sign(eta) * (0.4 * exp(-abs(eta) / 0.01) + 0.45 * exp(-abs(eta) / 0.1) + 0.35); % charlotte model
    
end