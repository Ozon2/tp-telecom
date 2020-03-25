% Troisième chaine étudiée : impact du choix de filtre de réception
close all
clear all

%% Initialisation des constantes
Fe = 12000;   % (Hz) fréquence d'échantillonnage
Te = 1/Fe;
Rs = 3000;    % (symboles) rythme symbole
alpha = 0.5;   % roll off du filtre de réception

%% Initialisation des constantes
Nb = 10000; % Nombre de bits
Ns = (1-alpha)/(1/Rs);    % Nombre d'échantillone par période symbole
SPAN = 1;
h = rcosdesign (alpha, SPAN, Ns, 'sqrt'); % Répertoire impulsionnelle du filtre de mise en forme
hr = fliplr(h);  %Filtre de réception adapté
t0 = 1;

%% Génération des bits et Mapping
bits = randi([0,1],1,Nb);
symboles = 3*(2*bits - 1);

peigne_dirac = kron(symboles, [1, zeros(1,Ns-1)]);
x = filter(h, 1, peigne_dirac);

%% Canal
r = x;

%% Reception
z = filter(hr, 1, r);

%% Echantillonnage
ze = z(t0:Ns:Ns*Nb);

%% Décision
bits_estimes = (ze > 0);
TEB = sum(bits(1:Nb) ~= bits_estimes)/Nb;
fprintf("Le TEB sans bruit vaut : %d \n", TEB);

%% Tracé de la chaine
figure;
plot((1:Nb*Ns)/Ns, peigne_dirac,'LineWidth',2);
hold on
plot((1:Nb*Ns)/Ns, x, 'r--','LineWidth',2);
plot((1:Nb*Ns)/Ns, z, 'b.-','LineWidth',2);
plot((t0:Ns:Nb*Ns)/Ns, ze ,'gp','LineWidth',3);
legend ("Peigne de Dirac", "x(t)", "z(t)", "ze(t)");
xlabel("Temps (en période T_s)");
title("Chaine de référence");
grid

%% Diagramme de l'oeil
eyediagram (z(length(h):Nb*Ns), 2*Ns, 2*Ns);
