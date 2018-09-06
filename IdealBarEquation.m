function [maxTotEnergy] = IdealBarEquation (fs, L, BC)
% fs = 44100;      % Sampling rate
% L = 1;          % String length
rho = 7850;     % Density of steel [kg/m^3]
k = 1 / fs;     % Time step
r = 0.0001;      % String radius
A = pi * r^2;   % Cross-sectional area
E = 2E11;       % Young's Modulus
inertia = pi / 4 * r^4; % Moment of inertia
kappa = (E * inertia) / (rho * A * L^4);  % Stiffness

h = sqrt(2 * kappa * k); % Grid spacing
muSq = (k * kappa / h^2)^2; % Courant number squared

P = 1/2; % plucking position

N = floor(L / h); % Number of grid-points

%% Raised cosine input
cosWidth = floor(N / 10);
raisedCos = 0.5 * (cos(pi:(2*pi)/cosWidth:3*pi) + 1);
PIdx = floor (P * N);

%% Initialise state vectors
u = zeros(N,1);
uPrev = zeros(N,1);
uNext = zeros(N, 1);

lengthSound = fs;

%% Matrix representation
if strcmp(BC, 'clamped')
    N = N - 4;
    Dxxxx = (sparse(3:N, 1:N-2, ones(1, N-2), N, N) + ...
            sparse(2:N, 1:N-1, -4 * ones(1, N-1), N, N) + ...
            sparse(1:N, 1:N, 6 * ones(1, N), N, N) + ...
            sparse(1:N-1, 2:N, -4 * ones(1, N-1), N, N) + ...
            sparse(1:N-2, 3:N, ones(1, N-2), N, N));
   
    I = sparse(1:N, 1:N, ones(1, N), N, N);
    N = N + 4;
    u(3:N-2) = rand(N-4,1) - 0.5; 
    uPrev(3:N-2) = rand(N-4,1) - 0.5;
end 
if strcmp(BC, 'ss')
        N = N - 2;
        Dxxxx = (sparse(3:N, 1:N-2, ones(1, N-2), N, N) + ...
            sparse(2:N, 1:N-1, -4 * ones(1, N-1), N, N) + ...
            sparse(1:N, 1:N, 6 * ones(1, N), N, N) + ...
            sparse(1:N-1, 2:N, -4 * ones(1, N-1), N, N) + ...
            sparse(1:N-2, 3:N, ones(1, N-2), N, N));
        Dxxxx(1, 1) = 5;
        Dxxxx(N, N) = 5;
        I = sparse(1:N, 1:N, ones(1, N), N, N);
        N = N + 2; 
        u(2:N-1) = rand(N-2,1) - 0.5; 
        uPrev(2:N-1) = rand(N-2,1) - 0.5;
end
if strcmp(BC, 'free')
    Dxx = (sparse(2:N, 1:N-1, [ones(1, N-2) 2], N, N) + ...
            sparse(1:N, 1:N, -2 * ones(1, N), N, N) + ...
            sparse(1:N-1, 2:N, [2 ones(1, N-2)], N, N));
%     Dxx2 = (sparse(2:N, 1:N-1, [ones(1, N-2) 1], N, N) + ...
%         sparse(1:N, 1:N, -1 * ones(1, N), N, N) + ...
%         sparse(1:N-1, 2:N, [1 ones(1, N-2)], N, N));
    Dxxxx = Dxx * Dxx;
%     Dxxxx2 = Dxx2 * Dxx2;
%     Dxxxx = (sparse(3:N, 1:N-2, ones(1, N-2), N, N) + ...
%                     sparse(2:N, 1:N-1, -4 * ones(1, N-1), N, N) + ...
%                     sparse(1:N, 1:N, 6 * ones(1, N), N, N) + ...
%                     sparse(1:N-1, 2:N, -4 * ones(1, N-1), N, N) + ...
%                     sparse(1:N-2, 3:N, ones(1, N-2), N, N));
%     Dxxxx(1, 1:3) = [6 -8 2];
%     Dxxxx(2, 1:4) = [-4 7 -4 1];
%     Dxxxx(N, N-2:N) = [2 -8 6];
%     Dxxxx(N-1, N-3:N) = [1 -4 7 -4];
    I = sparse(1:N, 1:N, ones(1, N), N, N);
%     u = rand(N,1) - 0.5; 
%     uPrev = rand(N,1) - 0.5;
    scalefac = ones(N,1);
    scalefac(1) = 0.5;
    scalefac(N) = 0.5;  
%     u(38:38+cosWidth) = raisedCos;
    u(1:N) = rand(N,1) - 0.5;
    uPrev = u;
    u2 = u;
    uPrev2 = u;
end

B = 2 * I - muSq * Dxxxx;

kinEnergy = zeros(lengthSound,1);
potEnergy = zeros(lengthSound,1);

if strcmp(BC, "clamped")
    for n = 1 : lengthSound
        uNext(3:N-2) = B * u(3:N-2) - uPrev(3:N-2);
        kinEnergy(n) = 1 / 2 * h * sum((1 / k * (u - uPrev)).^2);
        potEnergy(n) = kappa^2 / 2 * 1/h^3 * sum((u(3:N) - 2 * u(2:N-1) + u(1:N-2))...
            .* (uPrev(3:N) - 2 * uPrev(2:N-1) + uPrev(1:N-2)));
        uPrev = u;
        u = uNext;
%         plot(u)
%         drawnow;
    end
end

if strcmp(BC, "ss")
    for n = 1 : lengthSound
        uNext(2:N-1) = B * u(2:N-1) - uPrev(2:N-1);
        kinEnergy(n) = 1 / 2 * h * sum((1 / k * (u - uPrev)).^2);
        potEnergy(n) = kappa^2 / 2 * 1/h^3 * (sum((u(3:N) - 2 * u(2:N-1) + u(1:N-2))...
            .* (uPrev(3:N) - 2 * uPrev(2:N-1) + uPrev(1:N-2)))...
            + (- 2 * u(1)) * (- 2 * uPrev(1))...
            + (-2 * u(N)) * (-2 * uPrev(N)));
        uPrev = u;
        u = uNext;
    end
end

if strcmp(BC, "free")
    for n = 1 : lengthSound
        uNext = B * u - uPrev;
%         uNext2 = B2 * u2 - uPrev2;
        kinEnergy(n) = h / (2 * k^2) * sum(scalefac .* (u - uPrev).^2);
        potEnergy(n) = kappa^2 / (2 * h^3) * (sum((u(3:N) - 2 * u(2:N-1) + u(1:N-2))...
            .* (uPrev(3:N) - 2 * uPrev(2:N-1) + uPrev(1:N-2)))...
            + 0.5 * (-2 * u(1) + 2 * u(2)) * (-2 * uPrev(1) + 2 * uPrev(2))...
            + 0.5 * (-2 * u(N) + 2 * u(N-1)) * (-2 * uPrev(N) + 2 * uPrev(N-1)));
        if mod(n,100) == 0
%             clf
%         plot(kinEnergy(1:n)); hold on;
%         plot(potEnergy(1:n)); hold on;
% %         plot(kinEnergy(1:n) + potEnergy(1:n));
% %         plot(u); hold on;
%             plot(u);
%             drawnow;
        end
        uPrev = u;
        u = uNext;
%         uPrev2 = u2;
%         u2 = uNext2;
    end
end
totEnergy = kinEnergy + potEnergy;
totEnergy = (totEnergy-totEnergy(1))/totEnergy(1);
maxTotEnergy = max(totEnergy) - min(totEnergy);
plot(totEnergy);
