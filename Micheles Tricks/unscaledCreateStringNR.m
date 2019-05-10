function [B, C, N, h, Dxx, Dxxxx, s0, s1, bB, bC] = unscaledCreateStringNR(rho, A, T, E, I, L, s0, s1, k)
   
    kappa = sqrt(E * I / (rho * A));
    c = sqrt(T/(rho * A));

    h = sqrt((c^2 * k^2 + 4 * s1 * k + sqrt((c^2 * k^2 + 4 * s1 * k)^2 + 16 * kappa^2 * k^2)) / 2);
    N = floor(L/h); % Number of gridpoints
    h = L/N; % Recalculate gridspacing
    
    s0 = s0 * rho * A;
    s1 = s1 * rho * A;
    Ndec = N;
    N = N - 4;
    Ndec = Ndec - N;
    Dxxxx = (sparse(3:N, 1:N-2, ones(1, N-2), N, N) + ...
            sparse(2:N, 1:N-1, -4 * ones(1, N-1), N, N) + ...
            sparse(1:N, 1:N, 6 * ones(1, N), N, N) + ...
            sparse(1:N-1, 2:N, -4 * ones(1, N-1), N, N) + ...
            sparse(1:N-2, 3:N, ones(1, N-2), N, N));
    Dxx =   (sparse(2:N, 1:N-1, ones(1, N-1), N, N) + ...
            sparse(1:N, 1:N, -2 * ones(1, N), N, N) + ...
            sparse(1:N-1, 2:N, ones(1, N-1), N, N));
%     Dxxxx = Dxx*Dxx;
%     N = N - 2;
    B = (2 * rho * A / k^2 * eye(N) + T / h^2 * Dxx - E * I / h^4 * Dxxxx + 2 * s1 / (k * h^2) * Dxx) / (rho * A / k^2 + s0 / k);
    C = ((-rho * A / k^2 + s0 / k) * eye(N) - 2 * s1 / (k*h^2) * Dxx) / (rho * A / k^2 + s0/k);
    
    kOh = (kappa/h^2)^2;
    gOh = (c/h)^2;

    phi = (2 * s1)/(k * h^2);
    
    bB = (sparse(3:N, 1:N-2, kOh * ones(1, N-2), N, N) + ...
                sparse(2:N, 1:N-1, (-4 * kOh - gOh - phi) * ones(1, N-1), N, N) + ...
                sparse(1:N, 1:N, (6*kOh + 2 * gOh - 2/k^2 + 2 * phi) * ones(1, N), N, N) + ...
                sparse(1:N-1, 2:N, (-4 * kOh - gOh - phi) * ones(1, N-1), N, N) + ...
                sparse(1:N-2, 3:N, kOh * ones(1, N-2), N, N));
    bC = sparse(2:N, 1:N-1, phi * ones(1, N-1), N, N) + ...
                    sparse(1:N, 1:N, (2/k^2 - 2 * phi) * ones(1, N), N, N) + ...
                    sparse(1:N-1, 2:N, phi * ones(1, N-1), N, N);
    %     C = -((1 - s0 * k) * eye(N) + 2 * s1 * k * Dxx(2:end-1,2:end-1) / h^2) / (1 + s0 * k);
    
    N = N + 4;
    matVec = Ndec / 2 + 1 : N - Ndec / 2;
    Dxx2 = zeros(N);
    Dxx2(matVec, matVec) = Dxx;
    Dxx = Dxx2;
    
    Dxxxx2 = zeros(N);
    Dxxxx2(matVec, matVec) = Dxxxx;
    Dxxxx = Dxxxx2;
    
    B2 = zeros(N);
    B2(matVec, matVec) = B;
    B = B2;

    C2 = zeros(N);
    C2(matVec, matVec) = C;
    C = C2;
    
    bB2 = zeros(N);
    bB2(matVec, matVec) = bB;
    bB = bB2;

    bC2 = zeros(N);
    bC2(matVec, matVec) = bC;
    bC = bC2;
end