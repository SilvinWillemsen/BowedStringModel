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
cosWidth = floor(N / 5);
raisedCos = 0.5 * (cos(pi:(2*pi)/cosWidth:3*pi) + 1);
PIdx = floor (P * N);

%% Initialise state vectors
u = rand(N,1) - 0.5; 
uPrev = rand(N,1) - 0.5;
uNext = zeros(N, 1);

u2 = zeros(N, 1);
u2(ceil(N * P - cosWidth / 2 : ceil(N * P + cosWidth / 2))) = raisedCos;
uPrev2 = u2;
uNext2 = zeros(N, 1);

u3 = zeros(N, 1);
u3(ceil(N * P - cosWidth / 2 : ceil(N * P + cosWidth / 2))) = raisedCos;
uPrev3 = u3;
uNext3 = zeros(N, 1);

u4 = zeros(N, 1);
u4(ceil(N * P - cosWidth / 2 : ceil(N * P + cosWidth / 2))) = raisedCos;
uPrev4 = u4;
uNext4 = zeros(N, 1);

lengthSound = fs * 5;
drawBar = true;
matrix = true;

ssBounds = true;
freeBounds = true;

bounds = 'clamped';
%% Matrix representation
if matrix
    if strcmp(bounds, 'clamped')
        Dxx = (sparse(2:N, 1:N-1, ones(1, N-1), N, N) + ...
            sparse(1:N, 1:N, -2 * ones(1, N), N, N) + ...
            sparse(1:N-1, 2:N, ones(1, N-1), N, N));
        Dxxxx = Dxx * Dxx;
        range = 3 : N - 2;
    end 
    if strcmp(bounds, 'ss')
%         N = N - 2;
            Dxxxx = (sparse(3:N, 1:N-2, ones(1, N-2), N, N) + ...
                sparse(2:N, 1:N-1, -4 * ones(1, N-1), N, N) + ...
                sparse(1:N, 1:N, 6 * ones(1, N), N, N) + ...
                sparse(1:N-1, 2:N, -4 * ones(1, N-1), N, N) + ...
                sparse(1:N-2, 3:N, ones(1, N-2), N, N));
            Dxxxx(2, 2) = 5;
            Dxxxx(N-1, N-1) = 5;

%         Dxx2 = (sparse(2:N, 1:N-1, ones(1, N-1), N, N) + ...
%                     sparse(1:N, 1:N, -2 * ones(1, N), N, N) + ...
%                     sparse(1:N-1, 2:N, ones(1, N-1), N, N));
%         Dxxxx2 = Dxx2 * Dxx2;
%         N = N + 2;
%             Dxxxx([1 end],[1 end]) = 5;
        range = 2 : N - 1;    
    end
    if strcmp(bounds, 'free')
        Dxx = (sparse(2:N, 1:N-1, [ones(1, N-2) 2], N, N) + ...
                sparse(1:N, 1:N, -2 * ones(1, N), N, N) + ...
                sparse(1:N-1, 2:N, [2 ones(1, N-2)], N, N));
            Dxxxx = Dxx * Dxx;
            Dxxxx2 = (sparse(3:N, 1:N-2, ones(1, N-2), N, N) + ...
                sparse(2:N, 1:N-1, -4 * ones(1, N-1), N, N) + ...
                sparse(1:N, 1:N, 6 * ones(1, N), N, N) + ...
                sparse(1:N-1, 2:N, -4 * ones(1, N-1), N, N) + ...
                sparse(1:N-2, 3:N, ones(1, N-2), N, N));
            Dxxxx2(1, 1:3) = [6 -8 2];
            Dxxxx2(2, 1:4) = [-4 7 -4 1];
            Dxxxx2(N, N-2:N) = [2 -8 6];
            Dxxxx2(N-1, N-3:N) = [1 -4 7 -4];
        range = 1:N;
    end
    I = sparse(1:N, 1:N, ones(1, N), N, N);
    B = 2 * I - muSq * Dxxxx;
    out = zeros(lengthSound, 1);

    kinEnergy = zeros(lengthSound,1);
    potEnergy = zeros(lengthSound,1);
    for n = 1 : lengthSound
        uNext(range) = B(range, range) * u(range) - uPrev(range);
        kinEnergy(n) = 1 / 2 * h * sum((1 / k * (u - uPrev)).^2);
        potEnergy(n) = kappa^2 / 2 * 1/h^3 * sum((u(3:N) - 2 * u(2:N-1) + u(1:N-2))...
                .* (uPrev(3:N) - 2 * uPrev(2:N-1) + uPrev(1:N-2)));
        uPrev = u;
        u = uNext;
        if mod(n,1) == 0 && drawBar
            plot(uNext);
%             ylim([-1 1]);
            drawnow;
        end
    end
    totEnergy = kinEnergy + potEnergy;
    totEnergy = (totEnergy-totEnergy(1))/totEnergy(1);
    plot(totEnergy);
%% For-loop representation
else
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

                uNext3(l) = 2 * u3(l) - muSq * (6 * u3(l) - ...
                    4 * (u3(l + 1) + u3(l - 1)) + ...
                    (u3(l - 2) + u3(l + 2))) - uPrev3(l);

            else
                if l == 2 && ssBounds
                    uNext2(2) = 2 * u2(2) - muSq * (5 * u2(l) - 4 * u2(l + 1) + u2(l + 2)) - uPrev2(2);
                end

                if l == N - 1 && ssBounds
                    uNext2(l) = 2 * u2(l) - muSq * (5 * u2(l) - 4 * u2(l - 1) + u2(l - 2)) - uPrev2(l);
                end

                if l == 1 && freeBounds
                    uNext3(l) = 2 * u3(l) - muSq * (6 * u3(l) - 8 * u3(l + 1) + 2 * u3(l + 2)) - uPrev3(l);
                end
                if l == 2 && freeBounds
                    uNext3(l) = 2 * u3(l) - muSq * (7 * u3(l) - 4 * (u3(l - 1) + u3(l + 1)) + u3(l + 2)) - uPrev3(l);
                end
                if l == N && freeBounds
                    uNext3(l) = 2 * u3(l) - muSq * (6 * u3(l) - 8 * u3(l - 1) + 2 * u3(l - 2)) - uPrev3(l);
                end
                if l == N - 1 && freeBounds
                    uNext3(l) = 2 * u3(l) - muSq * (7  * u3(l) - 4 * (u3(l + 1) + u3(l - 1)) + u3(l - 2)) - uPrev3(l);
                end
            end
            if l > 1 && l < N
                potEnergy(n) = potEnergy(n) +  kappa^2 / 2 * 1/h^3 ...
                    * (u(l+1) - 2 * u(l) + u(l-1)) * (uPrev(l+1) - 2 * uPrev(l) + uPrev(l-1));
            end
            if l == 1
                potEnergy(n) = potEnergy(n) + kappa^2 / 2 * 1/h^3 ...
                    * (2 * u(l+1) - 2 * u(l)) * (2 * uPrev(l+1) - 2 * u(l));
            end
            if l == N
                potEnergy(n) = potEnergy(n) + kappa^2 / 2 * 1/h^3 ...
                    * (-2 * u(l) + 2 * u(l-1)) * (-2 * uPrev(l) + 2 * u(l-1));
            end
            kinEnergy(n) = kinEnergy(n) + 1 / 2 * h * ((1 / k * (u(l) - uPrev(l)))^2);
        end
%         uNext4(range) = B(range, range) * u4(range) - uPrev4(range);
%         kinEnergy2(n) = sum (1 / 2 * h * ((1 / k * (u4 - uPrev4)).^2));
%         potEnergy2(n) = sum(kappa^2 / 2 * 1/h^3 * (u4(3:N) - 2 * u4(2:N-1) + u4(1:N-2))...
%                     .* (uPrev4(3:N) - 2 * uPrev4(2:N-1) + uPrev4(1:N-2)));
%             totEnergy = kinEnergy + potEnergy;
%             totEnergy = (totEnergy-totEnergy(1))/totEnergy(1);
%             totEnergy2 = kinEnergy2 + potEnergy2;
%             totEnergy2 = (totEnergy2-totEnergy2(1))/totEnergy2(1);

%             plot(totEnergy(1:n));
%             hold on;
%             plot(totEnergy2(1:n));
%             drawnow;
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
        
        uPrev4 = u4;
        u4 = uNext4;

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
end