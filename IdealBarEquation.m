% clear all;
close all;
clc;

fs = 10000;      % Sampling rate
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

N = floor(L / h) + 1; % Number of grid-points

%% Raised cosine input
cosWidth = N / 20;
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

lengthSound = fs;
drawBar = false;
matrix = true;

%% Matrix representation
if matrix
    Dxx = (sparse(2:N, 1:N-1, ones(1, N-1), N, N) + ...
        sparse(1:N, 1:N, -2 * ones(1, N), N, N) + ...
        sparse(1:N-1, 2:N, ones(1, N-1), N, N));
    Dxxxx = Dxx * Dxx;

    I = sparse(1:N, 1:N, ones(1, N), N, N);
    B = 2 * I - muSq * Dxxxx;
    range = 3 : N - 2;
    out = zeros(lengthSound, 1);
    Buse = B(range, range);
    energyRange = 2:N - 1;
    kinEnergy = zeros(lengthSound,1);
    potEnergy = zeros(lengthSound,1);
    
    for t = 1 : lengthSound
        uNext(range) = Buse * u(range) - uPrev(range);
        out(t) = uNext(length(uNext) - PIdx);
        
        kinEnergy(t) = 1 / 2 * sum (h * ((1 / k * (u(energyRange) - uPrev(energyRange))).^2)); 
        potEnergy(t) = kappa^2 / (2 * h^3) * sum (...
            (u(energyRange + 1) - 2 * u(energyRange) + u(energyRange - 1)) .* ...
            (uPrev(energyRange + 1) - 2 * uPrev(energyRange) + uPrev(energyRange - 1)));
        uPrev = u;
        u = uNext;


        if mod(t,10) == 1 && drawBar
            clf
            plot(uNext);
            hold on;
            plot(uNext2);
            drawnow;

        end
    end
    totEnergy = kinEnergy + potEnergy;
    totEnergy = (totEnergy-totEnergy(1))/totEnergy(1);
    plot(totEnergy);
else
    %% For-loop representation

    kinEnergy = zeros(lengthSound,1);
    potEnergy = zeros(lengthSound,1);
    for n = 1 : lengthSound
        for l = 2 : length(u) - 1
            if l ~= 2 && l ~= length(u) - 1
                uNext(l) = (2 - 6 * muSq) * u(l) + ...
                            4 * muSq * (u(l + 1) + u(l - 1)) - ...
                            muSq * (u(l - 2) + u(l + 2)) - uPrev(l);
            end
            kinEnergy(n) = kinEnergy(n) + 1 / 2 * h * ((1 / k * (u(l) - uPrev(l)))^2);
            potEnergy(n) = potEnergy(n) +  kappa^2 / 2 * 1/h^3 ...
                * (u(l+1) - 2 * u(l) + u(l-1)) * (uPrev(l+1) - 2 * uPrev(l) + uPrev(l-1));
%                 2/h^3 * (u(2) - u(1)) * (uPrev(2) - uPrev(1))...
%                 + 2/h^3 * (-u(N) + u(N-1)) * (-uPrev(N) + uPrev(N-1))
    %         potEnergy2(n) = potEnergy2(n) +  kappa^2 / 2 * 1/h^3 ...
    %             * (u(l+1) - 2 * u(l) + u(l-1)) * (uPrev(l+1) - 2 * uPrev(l) + uPrev(l-1)) ...
    %             + 2/h^3 * (u(2) - u(1)) * (uPrev(2) - uPrev(1))...
    %             + 2/h^3 * (-u(N) + u(N-1)) * (-uPrev(N) + uPrev(N-1));
    %         


    %         potEnergy(n) = potEnergy(n) +  kappa^2 / 2 * h...
    %             * (1 / (2 * h^2) * (u(l + 1) - 2 * u(l) + u(l - 1) + uPrev(l + 1) - 2 * uPrev(l) + uPrev(l - 1)))^2 ...
    %             - (k^2 * kappa^2) / 8 * h...
    %             * h  * (1 / (h^2 * k) * (u(l + 1) - 2 * u(l) + u(l-1) - uPrev(l + 1) + 2 * uPrev(l) - uPrev(l - 1)))^2;
    %         potEnergy2(n) = potEnergy2(n) + ...
    %             kappa^2 / 2 * h...
    %             * (1 / (2 * h^2) * (u(l + 1) - 2 * u(l) + u(l - 1) + uPrev(l + 1) - 2 * uPrev(l) + uPrev(l - 1)))^2 ...
    %             - (k^2 * kappa^2) / 8 ...
    %             * h * (1 / (h^2 * k) * (u(l + 1) - uPrev(l + 1) - 2 * (u(l) - uPrev(l)) + u(l - 1) - uPrev(l - 1)))^2;

        end
        if mod(n,2) == 1 && drawBar
            plot(uNext);
            ylim([-1 1]);
            drawnow;
        end
        out(n) = uNext(PIdx);
        uPrev = u;
        u = uNext;
    end
    % totEnergy = kinEnergy + potEnergy;
    totEnergy = kinEnergy + potEnergy;
    % totEnergy = (totEnergy-totEnergy(1))/totEnergy(1);
    totEnergy = (totEnergy-totEnergy(1))/totEnergy(1);
    % plot(totEnergy);
    % hold on;
    plot(totEnergy);
%     hold on; plot(potEnergy2);
%     hold on; plot(kinEnergy2);
    % plot(out);
end