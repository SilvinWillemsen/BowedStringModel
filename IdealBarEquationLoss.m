clear all;
close all;
clc;

fs = 44100;      % Sampling rate
k = 1 / fs;     % Time step
L = 0.7;          % Bar length
rho = 7850;     % Density of steel [kg/m^3]
r = 0.008;      % String radius
A = pi * r^2;   % Cross-sectional area
E = 2e11;       % Young's Modulus
inertia = pi / 4 * r^4; % Moment of inertia
kappa = sqrt((E * inertia) / (rho * A * L^2));  % Stiffness

s0 = 0;
s1 = 0.0;
h = sqrt(2 * k * (s1^2 + sqrt(kappa^2 + s1^2)));
% muSq = (k * kappa / h^2)^2; % Courant number squared (always 0.25)
muSq = 0.25;
% lambdaSq = (k^2*kappa^2) / h^4;
lambdaSq = 0.25;
P = 1/2; % plucking position

N = floor(L / h); % Number of grid-points

%% Raised cosine input
cosWidth = floor(N / 2);
raisedCos = 0.5 * (cos(pi:(2*pi)/cosWidth:3*pi) + 1);
PIdx = floor (P * N);

%% Initialise state vectors
u = zeros(N, 1);
u(ceil(N * P - cosWidth / 2 : ceil(N * P + cosWidth / 2))) = raisedCos;
uPrev = u;
uNext = zeros(N, 1);

% u0 = u;
u0 = zeros(N,1);
uPrev0 = uPrev;
uNext0 = uNext;

u = zeros(N,1);
uPrev = zeros(N,1);
uNext = zeros(N,1);

u(1) = 0;
u(N) = 0;
uPrev = u;
% uPrev(1) = 0;
% uPrev(N) = 0;

u2 = zeros(N, 1);
u2(ceil(N * P - cosWidth / 2 : ceil(N * P + cosWidth / 2))) = raisedCos;
uPrev2 = u2;
uNext2 = zeros(N, 1);

u3 = zeros(N, 1);
u3(ceil(N * P - cosWidth / 2 : ceil(N * P + cosWidth / 2))) = raisedCos;
% u3 = rand(N, 1) - 0.5;
% u3 = rand(N, 1) - 0.5;
uPrev3 = u3;
uNext3 = zeros(N, 1);


u4 = zeros(N, 1);
u4(ceil(N * P - cosWidth / 2 : ceil(N * P + cosWidth / 2))) = raisedCos;
uPrev4 = u4;
uNext4 = zeros(N, 1);

lengthSound = fs;
drawBar = false;
matrix = false;

ssBounds = true;
freeBounds = true;

bounds = 'ss';
%% Matrix representation
if matrix
    if strcmp(bounds, 'clamped')
        N = N - 4;
        
        Dxxxx = (sparse(3:N, 1:N-2, ones(1, N-2), N, N) + ...
                sparse(2:N, 1:N-1, -4 * ones(1, N-1), N, N) + ...
                sparse(1:N, 1:N, 6 * ones(1, N), N, N) + ...
                sparse(1:N-1, 2:N, -4 * ones(1, N-1), N, N) + ...
                sparse(1:N-2, 3:N, ones(1, N-2), N, N));
        I = sparse(1:N, 1:N, ones(1, N), N, N);
        N = N + 4;
        u(37:37+cosWidth) = raisedCos;
        uPrev = u;
%         u(3:N-2) = rand(N-4,1) - 0.5; 
%         uPrev(3:N-2) = rand(N-4,1) - 0.5;
        
    elseif strcmp(bounds, 'ss')
        N = N - 2;
        Dxxxx = (sparse(3:N, 1:N-2, ones(1, N-2), N, N) + ...
            sparse(2:N, 1:N-1, -4 * ones(1, N-1), N, N) + ...
            sparse(1:N, 1:N, 6 * ones(1, N), N, N) + ...
            sparse(1:N-1, 2:N, -4 * ones(1, N-1), N, N) + ...
            sparse(1:N-2, 3:N, ones(1, N-2), N, N));
        Dxxxx(1, 1) = 5;
        Dxxxx(N, N) = 5;
        I = sparse(1:N, 1:N, ones(1, N), N, N);
%         Dxx2 = (sparse(2:N, 1:N-1, ones(1, N-1), N, N) + ...
%                     sparse(1:N, 1:N, -2 * ones(1, N), N, N) + ...
%                     sparse(1:N-1, 2:N, ones(1, N-1), N, N));
%         Dxxxx2 = Dxx2 * Dxx2;
        N = N + 2;
%             Dxxxx([1 end],[1 end]) = 5;
        u(2:N-1) = rand(N-2,1) - 0.5; 
        uPrev(2:N-1) = rand(N-2,1) - 0.5;   
    elseif strcmp(bounds, 'free')
        Dxx = (sparse(2:N, 1:N-1, [ones(1, N-2) 2], N, N) + ...
                sparse(1:N, 1:N, -2 * ones(1, N), N, N) + ...
                sparse(1:N-1, 2:N, [2 ones(1, N-2)], N, N));
%         Dxx = (sparse(2:N, 1:N-1, [ones(1, N-2) 2], N, N) + ...
%             sparse(1:N, 1:N, -2 * ones(1, N), N, N) + ...
%             sparse(1:N-1, 2:N, [2 ones(1, N-2)], N, N));
            Dxxxx = Dxx * Dxx;
%             Dxxxx = (sparse(3:N, 1:N-2, ones(1, N-2), N, N) + ...
%                 sparse(2:N, 1:N-1, -4 * ones(1, N-1), N, N) + ...
%                 sparse(1:N, 1:N, 6 * ones(1, N), N, N) + ...
%                 sparse(1:N-1, 2:N, -4 * ones(1, N-1), N, N) + ...
%                 sparse(1:N-2, 3:N, ones(1, N-2), N, N));
%             Dxxxx(1, 1:3) = [6 -8 2];
%             Dxxxx(2, 1:4) = [-4 7 -4 1];
%             Dxxxx(N, N-2:N) = [2 -8 6];
%             Dxxxx(N-1, N-3:N) = [1 -4 7 -4];
        I = sparse(1:N, 1:N, ones(1, N), N, N);
        u = rand(N,1) - 0.5; 
        uPrev = rand(N,1) - 0.5;
    end
    
    B = 2 * I - muSq * Dxxxx;
    out = zeros(lengthSound, 1);

    kinEnergy = zeros(lengthSound,1);
    potEnergy = zeros(lengthSound,1);
    if strcmp(bounds, 'clamped')
        for n = 1 : lengthSound
            uNext(3:N-2) = B * u(3:N-2) - uPrev(3:N-2);
            kinEnergy(n) = 1 / 2 * h * sum((1 / k * (u - uPrev)).^2);
            potEnergy(n) = kappa^2 / 2 * 1/h^3 * (sum((u(3:N) - 2 * u(2:N-1) + u(1:N-2))...
                    .* (uPrev(3:N) - 2 * uPrev(2:N-1) + uPrev(1:N-2))));%...
    %                 + (2 * u(2) - 2 * u(1)) * (2 * uPrev(2) - 2 * uPrev(1))...
    %                 + (-2 * u(N) + 2 * u(N-1)) * (-2 * uPrev(N) + 2 * uPrev(N-1)));
            out(n) = u(round(N/2));
            if mod(n,10) == 0
            plot(u)
            ylim([-1 1])
            drawnow;
            end
    uPrev = u;
            u = uNext;
        end
    elseif strcmp(bounds, 'ss')
        for n = 1 : lengthSound
            uNext(2:N-1) = B * u(2:N-1) - uPrev(2:N-1);
            kinEnergy(n) = 1 / 2 * h * sum((1 / k * (u - uPrev)).^2);
            potEnergy(n) = kappa^2 / 2 * 1/h^3 * (sum((u(3:N) - 2 * u(2:N-1) + u(1:N-2))...
                    .* (uPrev(3:N) - 2 * uPrev(2:N-1) + uPrev(1:N-2)))...
                    + (- 2 * u(1)) * (-2 * uPrev(1))...
                    + (-2 * u(N)) * (-2 * uPrev(N)));
            uPrev = u;
            u = uNext;
        end
    elseif strcmp(bounds, 'free')
        scalefac = ones(N,1);
        scalefac([1 N]) = 0.5;
        for n = 1 : lengthSound
            uNext = B * u - uPrev;
            kinEnergy(n) = 1 / 2 * h * sum(scalefac .* (1 / k * (u - uPrev)).^2);
            potEnergy(n) = kappa^2 / 2 * 1/h^3 * (sum((u(3:N) - 2 * u(2:N-1) + u(1:N-2))...
                    .* (uPrev(3:N) - 2 * uPrev(2:N-1) + uPrev(1:N-2)))...
                    + 0.5 * (- 2 * u(1) + 2 * u(2)) * (-2 * uPrev(1) + 2 * uPrev(2))...
                    + 0.5 * (-2 * u(N) + 2 * u(N-1)) * (-2 * uPrev(N) + 2 * uPrev(N-1)));
            uPrev = u;
            u = uNext;
        end
    end
    totEnergy = kinEnergy + potEnergy;
    totEnergy = (totEnergy-totEnergy(1))/totEnergy(1);
    plot(totEnergy);
%% For-loop representation
else
    kinEnergy = zeros(lengthSound,1);
    potEnergy = zeros(lengthSound,1);
%         u3 = rand(N,1) - 0.5; 
%         uPrev3 = rand(N,1) - 0.5;
    kinEnergy2 = zeros(lengthSound,1);
    potEnergy2 = zeros(lengthSound,1);
    
%     s1 = 0.00075;
   
    d = 1 + s0*k;
    T60 = 6 * log(10)/s0;
    for n = 1 : lengthSound
        for l = 1 : N
            if l > round(N/2) - 2 && l < round(N/2) + 2 && n < 11
                excitation = 1;
            else
                excitation = 0;
            end
              
            if l > 2 && l < N - 1
                uNext0(l) = (2 - 6 * muSq) * u0(l) + ...
                        4 * muSq * (u0(l + 1) + u0(l - 1)) - ...
                        muSq * (u0(l - 2) + u0(l + 2)) - uPrev0(l);
                uNext(l) = ((2 - 6 * muSq - 4 * s1 * (k / h^2)) * u(l) + ...
                    (4 * muSq + 2 * s1 * (k/h^2)) * (u(l+1) + u(l-1)) +...
                    (-1 + k * s0 + 4 * s1 * (k/h^2)) * uPrev(l) - ...
                    2 * s1 * (k/h^2) * (uPrev(l+1) + uPrev(l-1)) - ...
                    muSq * (u(l+2) + u(l-2))) / (1 + k*s0);
                
                uNext2(l) = ((2 - 6 * muSq - 4 * s1 * (k / h^2)) * u2(l) + ...
                    (4 * muSq + 2 * s1 * (k/h^2)) * (u2(l+1) + u2(l-1)) +...
                    (-1 + k * s0 + 4 * s1 * (k/h^2)) * uPrev2(l) - ...
                    2 * s1 * (k/h^2) * (uPrev2(l+1) + uPrev2(l-1)) - ...
                    muSq * (u2(l+2) + u2(l-2)) + excitation) / (1 + k*s0);
                
                uNext3(l) = ((2 - 6 * muSq - 4 * s1 * (k / h^2)) * u3(l) + ...
                    (4 * muSq + 2 * s1 * (k/h^2)) * (u3(l+1) + u3(l-1)) +...
                    (-1 + k * s0 + 4 * s1 * (k/h^2)) * uPrev3(l) - ...
                    2 * s1 * (k/h^2) * (uPrev3(l+1) + uPrev3(l-1)) - ...
                    muSq * (u3(l+2) + u3(l-2))) / (1 + k*s0);
                
            else
                if l == 2 && ssBounds
                    uNext2(l) = ((2 - 6 * muSq - 4 * s1 * (k / h^2)) * u2(l) + ...
                        (4 * muSq + 2 * s1 * (k/h^2)) * (u2(l+1) + u2(l-1)) +...
                        (-1 + k * s0 + 4 * s1 * (k/h^2)) * uPrev2(l) - ...
                        2 * s1 * (k/h^2) * (uPrev2(l+1) + uPrev2(l-1)) - ...
                        muSq * (u2(l+2) - u2(l)) + excitation) / (1 + k*s0);
                elseif l == N - 1 && ssBounds
                    uNext2(l) = ((2 - 6 * muSq - 4 * s1 * (k / h^2)) * u2(l) + ...
                        (4 * muSq + 2 * s1 * (k/h^2)) * (u2(l+1) + u2(l-1)) +...
                        (-1 + k * s0 + 4 * s1 * (k/h^2)) * uPrev2(l) - ...
                        2 * s1 * (k/h^2) * (uPrev2(l+1) + uPrev2(l-1)) - ...
                        muSq * (-u2(l) + u2(l-2)) + excitation) / (1 + k*s0);
                end
                
                if l == 1 && freeBounds
                    uNext3(l) = ((2 - 6 * muSq - 4 * s1 * (k / h^2)) * u3(l) + ...
                        (4 * muSq + 2 * s1 * (k/h^2)) * (2 * u3(l+1)) +...
                        (-1 + k * s0 + 4 * s1 * (k/h^2)) * uPrev3(l) - ...
                        2 * s1 * (k/h^2) * (2 * uPrev3(l+1)) - ...
                        muSq * (2 * u3(l+2))) / (1 + k*s0);
                end
                if l == 2 && freeBounds
                    uNext3(l) = ((2 - 6 * muSq - 4 * s1 * (k / h^2)) * u3(l) + ...
                        (4 * muSq + 2 * s1 * (k/h^2)) * (u3(l+1) + u3(l-1)) +...
                        (-1 + k * s0 + 4 * s1 * (k/h^2)) * uPrev3(l) - ...
                        2 * s1 * (k/h^2) * (uPrev3(l+1) + uPrev3(l-1)) - ...
                        muSq * (u3(l+2) + u3(l))) / (1 + k*s0);
                elseif l == N && freeBounds
                    uNext3(l) = ((2 - 6 * muSq - 4 * s1 * (k / h^2)) * u3(l) + ...
                        (4 * muSq + 2 * s1 * (k/h^2)) * (2 * u3(l-1)) +...
                        (-1 + k * s0 + 4 * s1 * (k/h^2)) * uPrev3(l) - ...
                        2 * s1 * (k/h^2) * (2 * uPrev3(l-1)) - ...
                        muSq * (2 * u3(l-2))) / (1 + k*s0);
                elseif l == N - 1 && freeBounds
                    uNext3(l) = ((2 - 6 * muSq - 4 * s1 * (k / h^2)) * u3(l) + ...
                        (4 * muSq + 2 * s1 * (k/h^2)) * (u3(l+1) + u3(l-1)) +...
                        (-1 + k * s0 + 4 * s1 * (k/h^2)) * uPrev3(l) - ...
                        2 * s1 * (k/h^2) * (uPrev3(l+1) + uPrev3(l-1)) - ...
                        muSq * (u3(l) + u3(l-2))) / (1 + k*s0);
                end
%                 if l == 1 && freeBounds
%                     uNext3(l) = 2 * u3(l) - muSq * (6 * u3(l) - 8 * u3(l + 1) + 2 * u3(l + 2)) - uPrev3(l);
%                 end
%                 if l == 2 && freeBounds
%                     uNext3(l) = 2 * u3(l) - muSq * (7 * u3(l) - 4 * (u3(l - 1) + u3(l + 1)) + u3(l + 2)) - uPrev3(l);
%                 end
%                 if l == N && freeBounds
%                     uNext3(l) = 2 * u3(l) - muSq * (6 * u3(l) - 8 * u3(l - 1) + 2 * u3(l - 2)) - uPrev3(l);
%                 end
%                 if l == N - 1 && freeBounds
%                     uNext3(l) = 2 * u3(l) - muSq * (7  * u3(l) - 4 * (u3(l + 1) + u3(l - 1)) + u3(l - 2)) - uPrev3(l);
%                 end

            end
            if strcmp(bounds, "clamped")
                if l > 1 && l < N
                    potEnergy(n) = potEnergy(n) +  kappa^2 / 2 * 1/h^3 ...
                        * (u(l+1) - 2 * u(l) + u(l-1)) * (uPrev(l+1) - 2 * uPrev(l) + uPrev(l-1));
                end
                kinEnergy(n) = kinEnergy(n) + 1 / 2 * h * ((1 / k * (u(l) - uPrev(l)))^2);
            end
            if strcmp(bounds, "ss")
                if l > 1 && l < N
                    potEnergy(n) = potEnergy(n) + kappa^2 / 2 * 1/h^3 ...
                        * (u2(l+1) - 2 * u2(l) + u2(l-1)) * (uPrev2(l+1) - 2 * uPrev2(l) + uPrev2(l-1));
                end
                if l == 1
                    potEnergy(n) = potEnergy(n) + kappa^2 / 2 * 1/h^3 ...
                        * (-2 * u2(1)) * (-2 * uPrev2(1));
                end
                if l == N
                    potEnergy(n) = potEnergy(n) + kappa^2 / 2 * 1/h^3 ...
                        * (-2 * u2(N)) * (-2 * uPrev2(N));
                end
                kinEnergy(n) = kinEnergy(n) + 1 / 2 * h * ((1 / k * (u2(l) - uPrev2(l)))^2);
  
            end
            
            if strcmp(bounds, "free")
                if l > 1 && l < N
                    potEnergy(n) = potEnergy(n) +  kappa^2 / 2 * 1/h^3 ...
                        * (u3(l+1) - 2 * u3(l) + u3(l-1)) * (uPrev3(l+1) - 2 * uPrev3(l) + uPrev3(l-1));
                    kinEnergy(n) = kinEnergy(n) + 1 / 2 * h * ((1 / k * (u3(l) - uPrev3(l)))^2);
                end
%                 if l == 2
%                     potEnergy(n) = potEnergy(n) + kappa^2 / 2 * 1/h^3 ...
%                         * 0.75 * (-2 * u3(1)+ 2 * u3(2)) * (-2 * uPrev3(1) + 2 * uPrev3(2));
%                     kinEnergy(n) = kinEnergy(n) + 1 / 2 * h * 0.75 * ((1 / k * (u3(l) - uPrev3(l)))^2);
%                 end
%                 if l == N-1
%                     potEnergy(n) = potEnergy(n) + kappa^2 / 2 * 1/h^3 ...
%                         * 0.75 * (-2 * u3(N) + 2 * u3(N-1)) * (-2 * uPrev3(N) + 2 * uPrev3(N-1));
%                     kinEnergy(n) = kinEnergy(n) + 1 / 2 * h * 0.75 * ((1 / k * (u3(l) - uPrev3(l)))^2);
%                 
%                 end
                if l == 1
                    potEnergy(n) = potEnergy(n) + kappa^2 / 2 * 1/h^3 ...
                        * 0.5 * (-2 * u3(1)+ 2 * u3(2)) * (-2 * uPrev3(1) + 2 * uPrev3(2));
                    kinEnergy(n) = kinEnergy(n) + 1 / 2 * h * 1 / 2 * ((1 / k * (u3(l) - uPrev3(l)))^2);
                end
                if l == N
                    potEnergy(n) = potEnergy(n) + kappa^2 / 2 * 1/h^3 ...
                        * 0.5 * (-2 * u3(N) + 2 * u3(N-1)) * (-2 * uPrev3(N) + 2 * uPrev3(N-1));
                    kinEnergy(n) = kinEnergy(n) + 1 / 2 * h * 1 / 2 * ((1 / k * (u3(l) - uPrev3(l)))^2);
                
                end
            end 
        end
%         clf
% if mod (n,3) == 0
%     clf
%     totEnergy = kinEnergy + potEnergy;
%     totEnergy = (totEnergy-totEnergy(1))/totEnergy(1);
%     plot(totEnergy(1:n)); 
%     hold on;
%     plot((potEnergy(1:n)-totEnergy(1))/totEnergy(1));
%     plot((kinEnergy(1:n)-totEnergy(1))/totEnergy(1));
%     drawnow;
% end
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
        if mod(n,100) == 0 && drawBar
%             clf;
%             plot(uNext); hold on;
%             plot(uNext2); hold on;
            plot(uNext3);
            ylim([-1 1]);
%             legend('Clamped', 'Simply Supported', 'Free');
            set(gca, 'FontSize', 15);
            drawnow;
        end
        out0(n) = uNext0(floor(length(uNext0) / 2));
        out(n) = uNext(floor(length(uNext) / 4));
        out2(n) = uNext2(floor(length(uNext) / 2));
        out3(n) = uNext3(floor(length(uNext) / 2));
%         clf;
%         if mod(n, 1) == 0
%                     plot(u);
% %             ylim([-1 1])
%             drawnow;
%         end
           
        uPrev0 = u0;
        u0 = uNext0;
        
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
%     totEnergy2 = kinEnergy2 + potEnergy2;
totEnergy = totEnergy(20:end);
    totEnergy = (totEnergy-totEnergy(1))/totEnergy(1);
%     totEnergy2 = (totEnergy2-totEnergy2(1))/totEnergy2(1);
    % plot(totEnergy);
    % hold on;
%     plot(totEnergy); 
%     hold on;
%     plot(totEnergy2);
    %     hold on; plot(potEnergy2);
    %     hold on; plot(kinEnergy2);
    plot(out3);
end