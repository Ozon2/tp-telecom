% Première chaine étudiée : chaine de référence
close all
clear all

%% Initialisation des constantes
Nb = 1000; % Nombre de bits
Ne = 8;    % Nombre d'échantillon par période symbole
h = ones(1, Ne); % Répertoire impulsionnelle du filtre de mise en forme
hr = fliplr(h);  %Filtre de réception adapté
t0 = length(h);

%% Génération des bits et Mapping
bits = randi([0,1],1,Nb);
symboles = 3*(2*bits - 1);

peigne_dirac = kron(symboles, [1, zeros(1,Ne-1)]);
x = filter(h, 1, peigne_dirac);

%% Densité spectrale de puissance
X = fft(x);
DSP = 1/(Nb*Ne) * abs(X).^2;
figure; semilogy(linspace(-0.5, 0.5, length(DSP)),fftshift(DSP));
title("DSP du signal transmis");
xlabel("Fréquence normalisée");
ylabel("DSP(f)");

%% Canal
r = x;

%% Reception
z = filter(hr, 1, r);

%% Echantillonnage
ze = z(t0:Ne:Ne*Nb);

%% Tracé de la chaine
figure;
plot((1:Nb*Ne)/Ne, peigne_dirac,'LineWidth',2);
hold on
plot((1:Nb*Ne)/Ne, x, 'r--','LineWidth',2);
plot((1:Nb*Ne)/Ne, z, 'b.-','LineWidth',2);
plot((t0:Ne:Nb*Ne)/Ne, ze ,'gp','LineWidth',3);
legend ("Peigne de Dirac", "x(t)", "z(t)", "ze(t)");
xlabel("Temps (en période T_s)");
title("Chaine de référence");
grid

%% Diagramme de l'oeil
eyediagram (z(length(h):Nb*Ne), 2*Ne, 2*Ne);

%% Décision
bits_estimes = (ze > 0);
TEB = sum(bits ~= bits_estimes)/Nb;
fprintf("Le TEB sans bruit vaut : %d \n", TEB);



%% TEB avec bruit
Eb_sur_N0_dB = linspace(0,6,100);
Eb_sur_N0 = 10.^(Eb_sur_N0_dB./10);
TEBs = zeros(1,length(Eb_sur_N0));
Pr = mean(abs(x).^2);
sigmas = Pr*Ne./(2*Eb_sur_N0);

TEB_theorique = 0;
p = 0.5; % probabilité d'une erreur (TODO)
Nelimit = Nb*p/(1-p); % Nombre d'erreurs à atteindre
for i = 1:length(sigmas)
    Nerr = 0; % Nombre d'erreurs
    nbEssais = 0;
    while (Nerr < Nelimit)S
        % Canal avec bruit AWGN
        r = x + sqrt(sigmas(i))*randn(1,length(x));
        % Réception
        z = filter(hr, 1, r);
        % Echantillonage
        ze = z(t0:Ne:Ne*Nb);
        % Décision
        bits_estimes = (ze > 0);
        NerrActuel =  sum(bits ~= bits_estimes);
        Nerr = Nerr + NerrActuel;
        TEBs(i) = TEBs(i) + NerrActuel/Nb;
        nbEssais = nbEssais + 1;
    end
    TEBs(i) = TEBs(i)/nbEssais;
end

figure;
semilogy(Eb_sur_N0_dB,TEBs);
title("TEB en fonction de (Eb/N0) (dB)");
xlabel("(Eb/N0) (dB)");
ylabel("TEB");