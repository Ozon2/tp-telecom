% Deuxième chaine étudiée : impact du choix de filtre de réception
close all
clear all

%% Initialisation des constantes
Nb = 10000; % Nombre de bits
Ns = 8;    % Nombre d'échantillon par période symbole
h = ones(1, Ns); % Répertoire impulsionnelle du filtre de mise en forme
hr1 = fliplr(h);  %Filtre de réception adapté chaine 1
hr2 = ones(1, Ns/2);  %Filtre de réception adapté chaine 2
t01 = length(h);
t02 = linspace(Ns/2, Ns, Nb);

%% Génération des bits et Mapping
bits = randi([0,1],1,Nb);
symboles = 3*(2*bits - 1);

peigne_dirac = kron(symboles, [1, zeros(1,Ns-1)]);
x = filter(h, 1, peigne_dirac);

%% Canal
r = x;

%% Reception
z = filter(hr2, 1, r);
z1 = filter(hr1, 1, r);


%% Densité spectrale de puissance
Z2 = fft(z);
DSP2 = 1/(Nb*Ns) * abs(Z2).^2;

Z1 = fft(z1);
DSP1 = 1/(Nb*Ns) * abs(Z1).^2;

figure; 
semilogy(linspace(-0.5, 0.5, length(DSP1)),fftshift(DSP1)); hold on;
semilogy(linspace(-0.5, 0.5, length(DSP2)),fftshift(DSP2));
title("DSP du signal transmis");
xlabel("Fréquence normalisée");
ylabel("DSP(f)");
legend("DSP chaine de référence","DSP chaine étudiée");

%% Echantillonnage
ze = z(t02:Ns:Ns*Nb);

%% Tracé de la chaine
figure;
plot((1:Nb*Ns)/Ns, peigne_dirac,'LineWidth',2);
hold on
plot((1:Nb*Ns)/Ns, x, 'r--','LineWidth',2);
plot((1:Nb*Ns)/Ns, z, 'b.-','LineWidth',2);
plot((t02:Ns:Nb*Ns)/Ns, ze ,'gp','LineWidth',3);
legend ("Peigne de Dirac", "x(t)", "z(t)", "ze(t)");
xlabel("Temps (en période T_s)");
title("Chaine de référence");
grid

%% Diagramme de l'oeil
eyediagram (z(length(h):Nb*Ns), 2*Ns, 2*Ns);

%% Décision
bits_estimes = (ze > 0);
TEB = sum(bits ~= bits_estimes)/Nb;
fprintf("Le TEB sans bruit vaut : %d \n", TEB);



%% TEB avec bruit
Eb_sur_N0_dB = linspace(0,6,100);
Eb_sur_N0 = 10.^(Eb_sur_N0_dB./10);
TEB1s = zeros(1,length(Eb_sur_N0));
TEB2s = zeros(1,length(Eb_sur_N0));
Pr = mean(abs(x).^2);
sigmas = Pr*Ns./(2*Eb_sur_N0);

Nelimite = 100;
for i = 1:length(sigmas)
    Nerr = 0;
    nbEssais = 0;
    while (Nerr < Nelimite)
        % Canal avec bruit AWGN
        r = x + sqrt(sigmas(i))*randn(1,length(x)); 
        % Réception
         z = filter(hr2, 1, r); 
        % Echantilonage
        ze = z(t02:Ns:Ns*Nb); 
        % Décision
        bits_estimes = (ze > 0);
        NerrActuel = sum(bits ~= bits_estimes);
        Nerr = Nerr + NerrActuel;
        nbEssais = nbEssais + 1;
    end
    TEB2s(i) = Nerr/(nbEssais*Nb);
end

for i = 1:length(sigmas)
    Nerr = 0;
    nbEssais = 0;
    while (Nerr < Nelimite)
        % Canal avec bruit AWGN
        r = x + sqrt(sigmas(i))*randn(1,length(x)); 
        % Réception
         z = filter(hr1, 1, r); 
        % Echantilonage
        ze = z(t01:Ns:Ns*Nb); 
        % Décision
        bits_estimes = (ze > 0);
        NerrActuel = sum(bits ~= bits_estimes);
        Nerr = Nerr + NerrActuel;
        nbEssais = nbEssais + 1;
    end
    TEB1s(i) = Nerr/(nbEssais*Nb);
end

figure;
semilogy(Eb_sur_N0_dB,TEB2s); hold on;
semilogy(Eb_sur_N0_dB,qfunc(sqrt(Eb_sur_N0)),'r+');
title("TEB en fonction de (Eb/N0) (dB)");
xlabel("(Eb/N0) (dB)");
ylabel("TEB");
legend("TEB simulé","TEB théorique");

figure;
semilogy(Eb_sur_N0_dB,TEB1s); hold on;
semilogy(Eb_sur_N0_dB,TEB2s);
title("TEB en fonction de (Eb/N0) (dB)");
xlabel("(Eb/N0) (dB)");
ylabel("TEB");
legend("TEB simulé chaine de référence","TEB simulé chaine 2");
