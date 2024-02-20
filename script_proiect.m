%% valori masurate
figure
plot(iesire_70)
hold on;
plot(iesire_130)
hold on; 
plot(iesire_145)
hold off;

%% identificare

%% comanda 70
s = tf('s');
Hs = (1.35/(1.42*s + 1)) * exp(-0.5*s);

figure
plot(iesire_70)
hold on
step(70*Hs)

%% comanda 130

Hs = (1.17/(1.80*s + 1)) * exp(-0.5*s);

figure
plot(iesire_130)
hold on
step(130*Hs)

%% comanda 145

Hs = (1.10/(1.74*s + 1)) * exp(-0.5*s);

figure
plot(iesire_145)
hold on
step(145*Hs)

%% Functia de transfer a procesului

Hs = (1.16/(2.1*s + 1)) * exp(-0.5*s);

figure
plot(iesire_130)
hold on; 
plot(iesire_145)

hold on
step(137.5*Hs)
hold off;

%% Regulator PI cu metoda predictor Smith
Hs = (1.16/(2.1*s + 1)) * exp(-0.5*s);

Hr = (0.86*(2.1*s+1))/(2.1*s);
Hr_disc=c2d(Hr,0.1);

% pentru verificare
Hd=(1.8*s+0.86)/(2.1*s-1+exp(-0.5*s));

figure
step(120*Hs)
hold on
step(120*(feedback(Hd*Hs,1)))
hold off

Hd_disc=tf(c2d(Hd,0.1));

%% Regulator RST
Hs = (1.16/(2.1*s + 1)) * exp(-0.5*s);
Hs_disc=c2d(Hs,0.1);

B = 0.05394;
A = [1 -0.9535];
d = 5;

omega=1.515;
teta=0.75;
H_impus = (omega*omega) / (omega*omega+2*omega*teta*s+s*s);
H_impusd = c2d(H_impus,0.1,'zoh');

P = [1 -1.778 0.7985];

M = zeros(7,7);
for i = 1:6
    M(i,i) = 1;
end

for i = 1:6
    M(i+1,i) = - 0.9535;
end

M(7,7) = 0.05394;

p = [P 0 0 0 0]'; 

pol = inv(M)*p;

S = [pol(1:6)]';
R = [pol(7:7)]';
T = H_impusd.den{1}/ sum(B);

tt = 4.6;
w = omega * tt;
H_impus = (w*w) / (w*w+2*w*teta*s+s*s);
H_impusd = c2d(H_impus,0.1,'zoh');

B1 = H_impusd.num{1};
A1 = H_impusd.den{1};

%% RST cu integrator fortat
Hs = (1.16/(2.1*s + 1)) * exp(-0.5*s);
Hs_disc=c2d(Hs,0.1);

B = 0.05394;
A = [1 -0.9535];
A_nou = conv(A,[1 -1]);
d = 5;

omega=1.515;
teta=0.75;
H_impus = (omega*omega) / (omega*omega+2*omega*teta*s+s*s);
H_impusd = c2d(H_impus,0.1,'zoh');

P = [1 -1.778 0.7985];

M = zeros(8,8);

for i = 1:7
    M(i,i) = 1;
end

for i = 1:6
    M(i+1,i) = -1.9535;
end

for i = 1:6
    M(i+2,i) = 0.9535;
end

M(7,7) = 0.0539;
M(8,8) = 0.0539;

p = [P 0 0 0 0 0]'; 

pol = inv(M)*p;

S1 = [pol(1:6)]';
S = conv(S1,[1 -1]);
R = [pol(7:8)]';
T = H_impusd.den{1}/ sum(B);

tt = 4.6;
w = omega * tt;

H_impus = (w*w) / (w*w+2*w*teta*s+s*s);
H_impusd = c2d(H_impus,0.1,'zoh');

B1 = H_impusd.num{1};
A1 = H_impusd.den{1};