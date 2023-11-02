% punkt pracy i stale obiektu
cw = 4200; % J/(kg*K)           cieplo wlasciwe
rho = 1000; % kg/m^3            gestosc
k = 160; % J/(s*m*K)            wspolczynnik przenikalnosci cieplnej
R = 0.56; % m                   promien zbiornika
d = 0.13; % m                   grubosc plaszcza
Vpl = 1.02; %m^3                objetosc plaszcza
Vflow0 = 0.005; %m^3/s          przeplyw w chwili zerowej
H0 = 1; %m                      wysokosc cieczy w zbiorniku
Pp = 1; % m^2                   pole podstawy zbiornika
Vzb = Pp * H0; % m^3            objetosc zbiotnika
T0 = 20; %st. C                 temperatura cieczy w zbiorniku
Tin0 = 70; %st. C               temperatura cieczy wplywajacej do plaszcza

% parametry transmitancji
% dTout(s)/dTin(s)
G11 = ((cw*rho*Vflow0-3.52*k*H0)*1)/(3.52*k*Vflow0);
T11 = Vpl/Vflow0;
T12 = (1*cw*rho)/(3.52*k);
% dT(s)/dTin(s)
G21 = 1;
T21 = (1*cw*rho)/(3.52*k);

% transmitancje
s = tf('s');
K1 = (s*G11+1)/((s*T11+1)*(s*T12+1));
K2 = G21/(T21*s+1);

% step response
figure
step(K1, 'r', 0:1:5000);
figure
step(K2, 'b', 0:1:500000);

% model simulink
sim('model_jedna_zmienna.mdl');
figure
plot(ans.tout, ans.simout);
hold on
grid on
legend('Tout(t)', 'T(t)');

% cw*rho*Vpl*(dTout(t)/dt) = cw*rho*Vflow(t)*Tin(t) -
% cw*rho*Vflow(t)*Tout(t) - 3.52*k*H0*(Tin(t) - T(t))

% 1*cw*rho*H0*(dT(t)/dt) = 3.52*k*H0*(Tin(t) - T(t))

Ts = 200000;                    % czas probkowania
dt = 0.1;                       % okres probkowania

Tout = zeros(1, Ts/dt);         % temperatura wyjsciowa
T = zeros(1, Ts/dt);            % temperatura w zbiorniku
Tin = 70*ones(1, Ts/dt);        % temperatura wejsciowa
Tin(50000:69999) = 50*ones(1, 20000);
Tin(70000:79999) = 60*ones(1, 10000);
Tin(120000:149999) = 50*ones(1, 30000);
Vflow = 0.0005*ones(1, Ts/dt);   % przeplyw
t = 1:1:(Ts/dt);                % wektor czasu

Tout_i = Tin(1);                % poczatkowa temperatura wyjsciowa
T_i = 20;                       % poczatkowa temperatura w zbiorniku

for i = 1:(Ts/dt)
    % w chwili i-tej
    T(i) = T_i;                 % temperatura w zbiorniku
    Tout(i) = Tout_i;           % temperatura wyjsciowa
    
    % obliczenie przyrostow temperatur
    dT = dt*(((3.52*k*H0)/(1*cw*rho*H0)) * (Tin(i) - T(i)));
    dTout = dt*((1/Vpl) * Vflow(i) * (Tin(i) - Tout(i))) - dT;
    
    % inkrementacja wartosci temperatur o ich przyrosty
    T_i = T_i + dT;
    Tout_i = Tout_i + dTout;
end

% wyswietlenie zaleznosci temperatur wyjsciowej i w zbiorniku od czasu
figure
plot(t, Tin, 'g');
hold on
plot(t, T, 'b');
plot(t, Tout, 'r');
xlabel('t (ms)');
ylabel('T (st. C)');
title('V* = 0.005m^3/s, k = 160J/s*m*K, H = 1m');
legend('Tin', 'T', 'Tout');
hold off

