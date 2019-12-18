clear all;
close all;

fs = 44100;
k = 1/fs;

lengthSound = fs;
init = true;
drawThings = true;
drawSpeed = 10;

%% Mass 1 variables
M1 = 0.001;
f1 = 500;
w1 = 2 * pi * f1;

%% Mass 2 variables
M2 = 0.01;
f2 = 200;
w2 = 2 * pi * f2;

%% Collision variables
K = 5e10;
alpha = 1.3;

%% Initialise state variables
u1Next = zeros(lengthSound, 1);
u1 = -0.999;
u1Prev = -1;

u2Next = zeros(lengthSound, 1);
u2 = 0;
u2Prev = 0;


%% NR tolerance
tol = 1e-13;

etaNext = u1 - u2;
for n = 1:lengthSound
    
    %% NR
    eps = 1;
    ii = 0;
    
    eta = u1 - u2;
    etaPrev = u1Prev - u2Prev;
    
    B = -2 * eta + etaPrev + k^2 * w1^2 * u1 - k^2 * w2^2 * u2;
    while eps > tol && ii < 100
        % variables for quotient rule
        f = subplus(etaNext)^(alpha+1) - subplus(etaPrev)^(alpha+1);
        df = (alpha+1) * subplus(etaNext)^alpha * 0.5 * (1+sign(etaNext));
        gg = etaNext - etaPrev;
        dgg = 1;

        G = etaNext + (K*k^2 / (M2*(alpha+1)) + K*k^2 / (M1*(alpha+1))) * (subplus(etaNext)^(alpha+1) - subplus(etaPrev)^(alpha+1)) / (etaNext - etaPrev) + B;
        Gdiff = 1 + (K*k^2 / (M2*(alpha+1)) + K*k^2 / (M1*(alpha+1))) * (df * gg - f * dgg) / gg^2;

        % NR for etaNext
        etaNextNRNew = etaNext - G / Gdiff;
        eps = abs(etaNextNRNew - etaNext);
        etaNext = etaNextNRNew;
        ii = ii + 1;
    end
    disp("NR iterations: " + ii)
    u1Next(n) = 2 * u1 - u1Prev - k^2 * w1^2 * u1 - K * k^2 / (M1 * (alpha+1)) * (subplus(etaNext)^(alpha+1) - subplus(etaPrev)^(alpha+1)) / (etaNext - etaPrev);
    u2Next(n) = 2 * u2 - u2Prev - k^2 * w2^2 * u2 + K * k^2 / (M2 * (alpha+1)) * (subplus(etaNext)^(alpha+1) - subplus(etaPrev)^(alpha+1)) / (etaNext - etaPrev);
    etaNext - (u1Next(n) - u2Next(n))
    % energy
    H1mass(n) = M1 / 2 * (1/k * (u1 - u1Prev))^2 + M1 * w1^2 / 2 * u1 * u1Prev;
    H2mass(n) = M2 / 2 * (1/k * (u2 - u2Prev))^2 + M2 * w2^2 / 2 * u2 * u2Prev;
    H1coll(n) = K / (alpha + 1) * 0.5 * (subplus(eta)^(alpha+1) + subplus(etaPrev)^(alpha+1));
    H1tot(n) = H1mass(n) + H2mass(n) + H1coll(n);
   
    %% Update states
    u1Prev = u1;
    u1 = u1Next(n);
    
    u2Prev = u2;
    u2 = u2Next(n);

    %% visualise stuff
    if drawThings && n > 10 && mod(n, drawSpeed) == 0
        subplot(2,1,1)
        hold off;
        plot(u1Next(1:n));
        hold on;
        plot(u2Next(1:n));
        title("State of the system")
        subplot(2,1,2)
        plot(H1tot(1:n) / H1tot(1) - 1);
%         plot(H1tot(1:n))
%         hold off
%         plot(H1mass(1:n))
%         hold on
%         plot(H2mass(1:n))
%         plot(H1coll(1:n))
        title("Normalised energy NR")
        drawnow;
    end
end