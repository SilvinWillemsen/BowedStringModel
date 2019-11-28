% Tromba test

% S. Bilbao
% 12/11/2019

clear all
close all
drawSpeed = 5;
SR = 44100;
T = 10;
rho = 7850;
r0 = 0.0003;
L = 1;
xM = 0.7;
Tf = 2;

amp = -0.1;
x0 = 0.5;
wid = 0.1;

M = 0.001;
b = -0.01;
alphaM = 1.3;
alphaB = 1.3;
KM = 5e14;
KB = 0e10;

A = pi*r0^2;
c = sqrt(T/(rho*A));
k = 1/SR;
h = c*k;
N = floor(L/h);
h = L/N;

u0 = zeros(N-1,1);
uMy1Next = zeros(N-1, 1);

for qq=1:N-1
    xcur = qq*h;
    dist = abs(xcur-x0);
    if(dist<=wid)
        u0(qq) = 0.5*amp*(1+cos(pi*dist/wid));
        uMy1Next(qq) = 0.5*amp*(1+cos(pi*dist/wid));
    end
end

u1 = u0;
u2 = u0;
uMy1 = uMy1Next;
uMy1Prev = uMy1Next;

w1 = -0.04;
w2 = -0.04;
psiM1 = 0;
psiB1 = 0;

e = ones(N-1,1);
Dxx = (1/h^2)*spdiags([e -2*e e],-1:1,N-1,N-1);

xMint = floor(xM*N/L);
J = sparse(N-1,1);
J(xMint) = 1;

Nf = floor(SR*Tf);

uMy2 = -0.04;
uMy2Prev = -0.04;

etaPrev = uMy2 - uMy1(xMint);
psiMyPrev = 0;
eta = uMy2Prev - uMy1Prev(xMint);
etaNext = uMy2 - uMy1(xMint);

vec = 2:N-2;
etaPrevTest= uMy2 - uMy1(xMint);
etaTest = uMy2Prev - uMy1Prev(xMint);
etaNextTest = uMy2 - uMy1(xMint);
psiMyPrevTest = 0;
for n=1:Nf
    
    etaB = b-w1;
    etaM = w1-J'*u1;
    etaB_plus = 0.5*(etaB+abs(etaB));
    etaM_plus = 0.5*(etaM+abs(etaM));
    
%     etaPrevTest = etaTest;
%     etaTest = etaNextTest;
%     
    etaPrev = eta;
    eta = etaNext;

%     eta = uMy2 - uMy1(xMint);
    gB = sqrt(KB*(alphaB+1)/2)*etaB_plus^((alphaB-1)/2);
    gM = sqrt(KM*(alphaM+1)/2)*etaM_plus^((alphaM-1)/2);
    gMy = sqrt(KM * (alphaM+1) / 2) * subplus(eta)^((alphaM - 1)/2);
    
    if gM ~= 0
        disp("wait");
    end
    A0 = [1 -k^2*gM^2/(4*M); 0 1+(k^2/(4*M))*gM^2+(k^2/(4*rho*A*h))*gM^2*(J'*J)];
        
    b0 = [-(w1-w2)/k+(k/(2*M))*gM*psiM1;   (w1-w2)/k-(k/(2*M))*gM*psiM1-J'*(u1-u2)/k-(T*k/(2*rho*A))*J'*(Dxx*u1)-(k/(2*rho*A*h))*(J'*J)*gM*psiM1];
    v = A0^(-1)*b0; % vector of \delta_{t dot}\eta <- comment added by silvin
    %\delta_{t\dot}\eta_M test
    test = ((w1-w2)/k-J'*(u1-u2)/k-(k/(2*M))*gM*psiM1-(T*k/(2*rho*A))*J'*(Dxx*u1)-(k/(2*rho*A*h))*gM*psiM1) ...
      / (1+(k^2/(4*M))*gM^2+(k^2/(4*rho*A*h))*gM^2);
%     psiB = psiB1+k*gB*v(1);
    psiM = psiM1+k*gM*v(2);

%     test - v(2)
    u = 2*u1-u2+c^2*k^2*Dxx*u1+(k^2/(rho*A))*(1/h)*J*gM*0.5*(psiM+psiM1);
    w = 2*w1-w2-(k^2/M)*gM*0.5*(psiM+psiM1);%+(k^2/M)*gB*0.5*(psiB+psiB1);
    
%     
    Amat = [rho*A/k^2+gMy^2/(4*h),   -gMy^2/(4*h);...
                  -gMy^2/4,          M/k^2+gMy^2/4];
              
    vMy = [rho*A/k^2 * (2 * uMy1(xMint) - uMy1Prev(xMint)) + T / h^2 * (uMy1(xMint+1) - 2 * uMy1(xMint) + uMy1(xMint-1)) ...
            - (gMy^2 / 4 * etaPrev - psiMyPrev * gMy) / h;...
            M / k^2 * (2 * uMy2 - uMy2Prev) + gMy^2/4 * etaPrev - psiMyPrev * gMy];
%     
    solut = Amat\vMy;
    
    test2 = ((uMy2-uMy2Prev)/k-J'*(uMy1-uMy1Prev)/k-(k/(2*M))*gM*psiMyPrev-(T*k/(2*rho*A))*J'*(Dxx*uMy1)-(k/(2*rho*A*h))*gMy*psiMyPrev) ...
      / (1+(k^2/(4*M))*gMy^2+(k^2/(4*rho*A*h))*gMy^2);
    etaNext = solut(2)-solut(1);
    psiMy = psiMyPrev+k*gMy*test2;
%     psiMy = gMy / 2 * (etaNext - etaPrev) + psiMyPrev;
    uMy1Next = 2 * uMy1 - uMy1Prev + T * k^2 / (rho * A) * Dxx * uMy1 + (k^2/(rho*A*h)) * J * gMy * 0.5 * (psiMy + psiMyPrev);
    uMy2Next = 2 * uMy2 - uMy2Prev - (k^2 / M) * gMy * 0.5 * (psiMy + psiMyPrev);
%     uMy1Next = 2 * uMy1 - uMy1Prev + T * k^2 / (rho * A) * Dxx * uMy1 + (k^2/(rho*A*h)) * J * (gMy^2 / 4 * (etaNext - etaPrev) + psiMyPrev * gMy);
%     uMy2Next = 2 * uMy2 - uMy2Prev - (k^2 / M) * (gMy^2 / 4 * (etaNext - etaPrev) + psiMyPrev * gMy);
    diff(n) = uMy1Next(xMint) - solut(1);
    etaNext = uMy2Next - uMy1Next(xMint);
%     (gMy / 2 * (psiMy + psiMyPrev)) - (gMy^2 / 4 * (etaNext - etaPrev) + psiMyPrev * gMy) 
    
    HS(n) = 0.5*h*rho*A*sum((u-u1).^2)/k^2+0.5*T*sum((u(2:end)-u(1:end-1)).*(u1(2:end)-u1(1:end-1)))/h;
    HM(n) = 0.5*M*(w-w1)^2/k^2;
    HC(n) = 0.5*psiM^2;%+0.5*psiB^2;
    H(n) = 0.5*M*(w-w1)^2/k^2+0.5*h*rho*A*sum((u-u1).^2)/k^2+0.5*T*sum(([u;0]-[0;u]).*([u1;0]-[0;u1]))/h+0.5*psiM^2;%+0.5*psiB^2;
    HMy(n) = 0.5*M*(uMy2Next-uMy2)^2/k^2+0.5*h*rho*A*sum((uMy1Next-uMy1).^2)/k^2+0.5*T*sum(([uMy1Next;0]-[0;uMy1Next]).*([uMy1;0]-[0;uMy1]))/h+0.5*psiMy^2;
       
    if mod(n, drawSpeed) == 0
        subplot(5,1,1)
        hold off
        plot([1:N-1]'*h, u, 'k');
        hold on
        plot([1:N-1]'*h, uMy1Next, 'b');
        plot(xM*L, w,'r.', 'MarkerSize', 24);
        plot(xM*L, b, 'b.', 'MarkerSize', 24);
        subplot(5,1,2)
        if n > 10
            plot(H(10:n) / H(10) - 1)
        end
        subplot(5,1,3)
        hold off
        plot([1:N-1]'*h, uMy1, 'k');
        hold on
        plot(xM*L, uMy2,'r.', 'MarkerSize', 24);
        plot(xM*L, b, 'b.', 'MarkerSize', 24);
        subplot(5,1,4)
        if n > 10
            plot(HMy(10:n) / HMy(10) - 1)
        end
        subplot(5,1,5)
        plot(u - uMy1Next)
        drawnow
    end
    u2 = u1;
    u1 = u;
    
    w2 = w1;
    w1 = w;
    
%     psiB1 = psiB;
    psiM1 = psiM;
    
    uMy1Prev = uMy1;
    uMy1 = uMy1Next;
    
    uMy2Prev = uMy2;
    uMy2 = uMy2Next;
    
    psiMyPrev = psiMy;
end
