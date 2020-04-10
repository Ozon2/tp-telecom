% Première chaine étudiée : chaine de référence
close all
clear all

%% Initialisation des constantes
Nb = 10000; % Nombre de bits
Ns = 8;    % Nombre d'échantillon par période symbole
h = ones(1, Ns); % Répertoire impulsionNslle du filtre de mise en forme
hr = fliplr(h);  %Filtre de réception adapté
t0 = length(h);

%% Génération des bits et Mapping
bits = randi([0,1],1,Nb);
symboles = 2*bits - 1;

peigne_dirac = kron(symboles, [1, zeros(1,Ns-1)]);
x = filter(h, 1, peigne_dirac);

%% Canal
r = x;

%% Reception
z = filter(hr, 1, r);



%% Densité spectrale de puissance
Z = fft(z);
DSP = 1/(Nb*Ns) * abs(Z).^2;
figure; semilogy(linspace(-0.5, 0.5, length(DSP)),fftshift(DSP));
title("DSP du signal transmis");
xlabel("Fréquence normalisée");
ylabel("DSP(f)");

%% Echantillonnage
ze = z(t0:Ns:Ns*Nb);

%% Tracé de la chaine
figure;
plot((1:Nb*Ns)/Ns, peigne_dirac,'LineWidth',2);
hold on
plot((1:Nb*Ns)/Ns, x, 'r--','LineWidth',2);
plot((1:Nb*Ns)/Ns, z, 'b.-','LineWidth',2);
plot((t0:Ns:Nb*Ns)/Ns, ze ,'gp','LineWidth',3);
legend ("Peigne de Dirac", "x(t)", "z(t)", "ze(t)");
xlabel("Temps (en période T_s)");
title("chaine de référence");
grid

%% Diagramme de l'oeil
eyediagram (z(length(h):Nb*Ns)/24, 2*Ns, 2*Ns);

%% Décision
bits_estimes = (ze > 0);
TEB = sum(bits ~= bits_estimes)/Nb;
fprintf("Le TEB sans bruit vaut : %d \n", TEB);

%% TEB avec bruit
Eb_sur_N0_dB = linspace(0,6,50);
Eb_sur_N0 = 10.^(Eb_sur_N0_dB./10);
TEBs = zeros(1,length(Eb_sur_N0));
Pr = mean(abs(x).^2);
sigmas = Pr*Ns./(2*Eb_sur_N0);

Nelimite = 1000;
for i = 1:length(sigmas)
    Nerr = 0;
    nbEssais = 0;
    while (Nerr < Nelimite)
        % Canal avec bruit AWGN
        r = x + sqrt(sigmas(i))*randn(1,length(x));
        % Réception
         z = filter(hr, 1, r);
        % Echantilonage
        ze = z(t0:Ns:Ns*Nb);
        % Décision
        bits_estimes = (ze > 0);
        NerrActuel = sum(bits ~= bits_estimes);
        Nerr = Nerr + NerrActuel;
        nbEssais = nbEssais + 1;
    end
    TEBs(i) = Nerr/(nbEssais*Nb);
end

figure;
semilogy(Eb_sur_N0_dB,TEBs, 'r+'); hold on;
semilogy(Eb_sur_N0_dB,qfunc(sqrt(2*Eb_sur_N0)),'g');
title("TEB en fonction de (Eb/N0) (dB)");
legend("TEB expérimental","TEB théorique");
xlabel("(Eb/N0) (dB)");
ylabel("TEB");
legend("TEB simulé","TEB théorique");
