clear all;
close all;
clc;

fs = 44100;      % Sampling rate
k = 1 / fs;     % Time step
L = 1;          % String length
rho = 7850;     % Density of steel [kg/m^3]
r = 0.001;      % String radius
A = pi * r^2;   % Cross-sectional area
E = 2E11;       % Young's Modulus
inertia = pi / 4 * r^4; % Moment of inertia
kappa = sqrt((E * inertia) / (rho * A * L^4));  % Stiffness

h = sqrt(2 * kappa * k); % Grid spacing
muSq = (k * kappa / h^2)^2; % Courant number squared

P = 1/2; % plucking position

N = floor(L / h); % Number of grid-points

%% Raised cosine input
cosWidth = floor(N / 20);
raisedCos = 0.5 * (cos(pi:(2*pi)/cosWidth:3*pi) + 1);
PIdx = floor (P * N);

%% Initialise state vectors
u = zeros(N, 1);
u(floor(N * P - cosWidth / 2 : ceil(N * P + cosWidth / 2))) = raisedCos;
uPrev = u;
uNext = zeros(N, 1);

u2 = zeros(N, 1);
u2(floor(N * P - cosWidth / 2 : ceil(N * P + cosWidth / 2))) = raisedCos;
uPrev2 = u2;
uNext2 = zeros(N, 1);

u3 = zeros(N, 1);
u3(floor(N * P - cosWidth / 2 : ceil(N * P + cosWidth / 2))) = raisedCos;
uPrev3 = u3;
uNext3 = zeros(N, 1);

lengthSound = fs * 5;
drawBar = false;
matrix = true;

ssBounds = true;
freeBounds = false;

bounds = 'SS';
%% Matrix representation
if strcmp(bounds, 'clamped')
    Dxx = (sparse(2:N, 1:N-1, ones(1, N-1), N, N) + ...
        sparse(1:N, 1:N, -2 * ones(1, N), N, N) + ...
        sparse(1:N-1, 2:N, ones(1, N-1), N, N));
    Dxxxx = Dxx * Dxx;
    range = 3 : N - 2;
else 
    if strcmp(bounds, 'SS')
        N = N - 2;
            Dxxxx = (sparse(3:N, 1:N-2, ones(1, N-2), N, N) + ...
                sparse(2:N, 1:N-1, -4 * ones(1, N-1), N, N) + ...
                sparse(1:N, 1:N, 6 * ones(1, N), N, N) + ...
                sparse(1:N-1, 2:N, -4 * ones(1, N-1), N, N) + ...
                sparse(1:N-2, 3:N, ones(1, N-2), N, N));
            Dxxxx(1, 1) = 5;
            Dxxxx(N, N) = 5;

%         Dxx2 = (sparse(2:N, 1:N-1, ones(1, N-1), N, N) + ...
%                     sparse(1:N, 1:N, -2 * ones(1, N), N, N) + ...
%                     sparse(1:N-1, 2:N, ones(1, N-1), N, N));
%         Dxxxx2 = Dxx2 * Dxx2;
        N = N + 2;
%             Dxxxx([1 end],[1 end]) = 5;
        range = 2 : N - 1;    
    end
end
I = sparse(1:N, 1:N, ones(1, N), N, N);
B = 2 * I(range, range) - muSq * Dxxxx;
out = zeros(lengthSound, 1);

kinEnergy = zeros(lengthSound,1);
potEnergy = zeros(lengthSound,1);

if matrix
    for n = 1 : lengthSound
        uNext(range) = B * u(range) - uPrev(range);
        kinEnergy(n) = 1 / 2 * h * sum((1 / k * (u - uPrev)).^2);
        potEnergy(n) = kappa^2 / 2 * 1/h^3 * sum((u(3:N) - 2 * u(2:N-1) + u(1:N-2))...
                .* (uPrev(3:N) - 2 * uPrev(2:N-1) + uPrev(1:N-2)));
        uPrev = u;
        u = uNext;
        if mod(n,1) == 0 && drawBar
            plot(uNext);
            ylim([-1 1]);
            drawnow;
        end
    end
    totEnergy = kinEnergy + potEnergy;
    totEnergy = (totEnergy-totEnergy(1))/totEnergy(1);
    plot(totEnergy);
end
%% For-loop representation

kinEnergy = zeros(lengthSound,1);
potEnergy = zeros(lengthSound,1);

kinEnergy2 = zeros(lengthSound,1);
potEnergy2 = zeros(lengthSound,1);
for n = 1 : lengthSound
    for l = 1 : N
        if l > 2 && l < length(u) - 1
            uNext(l) = (2 - 6 * muSq) * u(l) + ...
                        4 * muSq * (u(l + 1) + u(l - 1)) - ...
                        muSq * (u(l - 2) + u(l + 2)) - uPrev(l);

            uNext2(l) = 2 * u2(l) - muSq * (6 * u2(l) - ...
                4 * (u2(l + 1) + u2(l - 1)) + ...
                (u2(l - 2) + u2(l + 2))) - uPrev2(l);

            uNext3(l) = (2 - 6 * muSq) * u3(l) + ...
                4 * muSq * (u3(l + 1) + u3(l - 1)) - ...
                muSq * (u3(l - 2) + u3(l + 2)) - uPrev3(l);

        else
            if l == 2 && ssBounds
                uNext2(2) = 2 * u2(2) - muSq * (5 * u2(l) - 4 * u2(l + 1) + u2(l + 2)) - uPrev2(2);
            end

            if l == N - 1 && ssBounds
                uNext2(l) = 2 * u2(l) - muSq * (5 * u2(l) - 4 * u2(l - 1) + u2(l - 2)) - uPrev2(l);
            end

            if l == 1 && freeBounds
                uNext3(1) = (2 - 6 * muSq) * u3(1) + 8 * muSq * u3(2) - 2 * muSq * u3(3) - uPrev3(1);
            end
            if l == 2 && freeBounds
                uNext3(2) = (2 - 7 * muSq) * u3(2) + 4 * muSq * (u3(1) + u3(3)) - muSq * u3(4) - uPrev3(2);
            end
            if l == N && freeBounds
                uNext3(N) = (2 - 6 * muSq) * u3(N) + 8 * muSq * u3(N - 1) - 2 * muSq * u3(N - 2) - uPrev3(N);
            end
            if l == N - 1 && freeBounds
                uNext3(N - 1) = (2 - 7 * muSq) * u3(N - 1) + 4 * muSq * (u3(N) + u3(N - 2)) - muSq * u3(N - 3) - uPrev3(N - 1);
            end
        end
        if l > 1 && l < length(u)
            potEnergy(n) = potEnergy(n) +  kappa^2 / 2 * 1/h^3 ...
                * (u(l+1) - 2 * u(l) + u(l-1)) * (uPrev(l+1) - 2 * uPrev(l) + uPrev(l-1));
        end
        kinEnergy(n) = kinEnergy(n) + 1 / 2 * h * ((1 / k * (u(l) - uPrev(l)))^2);
    end
    kinEnergy2(n) = sum (1 / 2 * h * ((1 / k * (u - uPrev)).^2));
    potEnergy2(n) = sum(kappa^2 / 2 * 1/h^3 * (u(3:N) - 2 * u(2:N-1) + u(1:N-2))...
                .* (uPrev(3:N) - 2 * uPrev(2:N-1) + uPrev(1:N-2)));
%         totEnergy = kinEnergy + potEnergy;
%         totEnergy = (totEnergy-totEnergy(1))/totEnergy(1);
%         totEnergy2 = kinEnergy2 + potEnergy2;
%         totEnergy2 = (totEnergy2-totEnergy2(1))/totEnergy2(1);

%         plot(potEnergy(1:n) - potEnergy2(1:n));
%         hold on;
%         plot(potEnergy2(1:n));
%         drawnow;
    if mod(n,1) == 0 && drawBar
        clf;
        plot(uNext); hold on;
        plot(uNext2); hold on;
        plot(uNext3);
        ylim([-1 1]);
        legend('Clamped', 'Simply Supported', 'Free');
        set(gca, 'FontSize', 15);
        drawnow;
    end
    out(n) = uNext(floor(length(uNext) / 2));
    out2(n) = uNext2(floor(length(uNext) / 2));
    out3(n) = uNext3(floor(length(uNext) / 2));
    uPrev = u;
    u = uNext;

    uPrev2 = u2;
    u2 = uNext2;

    uPrev3 = u3;
    u3 = uNext3;
    
end
totEnergy = kinEnergy + potEnergy;
totEnergy2 = kinEnergy2 + potEnergy2;
totEnergy = (totEnergy-totEnergy(1))/totEnergy(1);
totEnergy2 = (totEnergy2-totEnergy2(1))/totEnergy2(1);
% plot(totEnergy);
% hold on;
plot(totEnergy); hold on;
plot(totEnergy2);
%     hold on; plot(potEnergy2);
%     hold on; plot(kinEnergy2);
% plot(out);