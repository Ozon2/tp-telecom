% Première chaine étudiée : chaine de référence
close all
clear all

%% Initialisation des coNetantes
Nb = 100; % Nombre de bits
Ne = 10;    % Nombre d'échantillone par période symbole
SPAN = 3;
h = ones(1, Ne); % Répertoire impulsionnelle du filtre de mise en forme
% h = rcosdesign (0.5, SPAN, Ne, 'sqrt');
hr = fliplr(h);  %Filtre de réception adapté
t0 = length(h);

%% Génération des bits et Mapping
bits = randi([0,1],1,Nb);
symboles = 3*(2*bits - 1);

peigne_dirac = kron(symboles, [1, zeros(1,Ne-1)]);
x = filter(h, 1, peigne_dirac);

%% Canal
r = x;

%% Reception
z = filter(hr, 1, r);

%% Echantillonnage
ze = z(t0:Ne:Ne*Nb);

%% Décision
bits_estimes = (ze > 0);
TEB = sum(bits ~= bits_estimes)/Nb;
% TEB = sum(bits(1:Nb-SPAN) ~= bits_estimes)/Nb;

%% Tracé de la chaiNe
figure;
plot((1:Nb*Ne)/Ne, peigne_dirac,'LineWidth',2);
hold on
plot((1:Nb*Ne)/Ne, x, 'r--','LineWidth',2);
plot((1:Nb*Ne)/Ne, z, 'b.-','LineWidth',2);
plot((t0:Ne:Nb*Ne)/Ne, ze ,'gp','LineWidth',3);
legend ("Peigne de Dirac", "x(t)", "z(t)", "ze(t)");
xlabel("Temps (en période T_s)");
title("chaines de référence");
grid

%% Diagramme de l'oeil
eyediagram (z(length(h):Nb*Ne), 2*Ne, 2*Ne);














