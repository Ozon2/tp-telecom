% Première chaine étudiée : chaine de référence
close all
clear all

%% Première chaine
% Initialisation des constantes chaine 1
Nb = 10000; % Nombre de bits
Ns = 8;    % Nombre d'échantillon par période symbole
h = ones(1, Ns); % Répertoire impulsionNslle du filtre de mise en forme
hr = fliplr(h);  %Filtre de réception adapté
t0 = length(h);

% Génération des bits et Mapping chaine 1
bits = randi([0,1],1,Nb);
symboles = 2*bits - 1;
peigne_dirac = kron(symboles, [1, zeros(1,Ns-1)]);

x = filter(h, 1, peigne_dirac);

% Canal chaine 1
r = x;

% Reception chaine 1
z = filter(hr, 1, r);


% Densité spectrale de puissance chaine 1
Z = fft(z);
DSP1 = 1/(Nb*Ns) * abs(Z).^2;

% TEB avec bruit chaine 1
Eb_sur_N0_dB = linspace(0,6,50);
Eb_sur_N0 = 10.^(Eb_sur_N0_dB./10);
TEB1s = zeros(1,length(Eb_sur_N0));
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
    TEB1s(i) = Nerr/(nbEssais*Nb);
end


% Troisième chaine étudiée : impact du choix de filtre de réception
%% Initialisation des constantes 
Fe = 12000;   % (Hz) fréquence d'échantillonnage
Te = 1/Fe;
Rs = 3000;    % (symboles) rythme symbole
alpha = 0.5;   % roll off du filtre de réception

%% Initialisation des constantes
Nb = 1000; % Nombre de bits
Ns = Fe/Rs;    % Nombre d'échantillone par période symbole
SPAN = 3;
h = rcosdesign (alpha, SPAN, Ns, 'sqrt'); % Répertoire impulsionnelle du filtre de mise en forme
hr = fliplr(h);  %Filtre de réception adapté
t0 = SPAN*Ns+1;

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
figure; semilogy(linspace(-0.5, 0.5, length(DSP)),fftshift(DSP)); hold on;
semilogy(linspace(-0.5, 0.5, length(DSP1)),fftshift(DSP1));
title("DSP du signal transmis");
xlabel("Fréquence normalisée");
ylabel("DSP(f)");
legend("DSP chaine de référence","DSP chaine étudiée");

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
title("Chaine de référence");
grid

%% Diagramme de l'oeil
eyediagram (z(length(h):Nb*Ns), 2*Ns, 2*Ns);

%% Décision
bits_estimes = (ze > 0);
TEB = sum(bits(1:Nb-SPAN) ~= bits_estimes)/Nb;
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
        NerrActuel = sum(bits(1:Nb-SPAN) ~= bits_estimes);
        Nerr = Nerr + NerrActuel;
        nbEssais = nbEssais + 1;
    end
    TEBs(i) = Nerr/(nbEssais*Nb);
end

figure;
semilogy(Eb_sur_N0_dB,TEBs, 'r+'); hold on;
semilogy(Eb_sur_N0_dB,qfunc(sqrt(2*Eb_sur_N0)), 'g');
title("TEB en fonction de (Eb/N0) (dB)");
xlabel("(Eb/N0) (dB)");
ylabel("TEB");
legend("TEB simulé","TEB théorique");

figure;
semilogy(Eb_sur_N0_dB,TEBs); hold on;
semilogy(Eb_sur_N0_dB,TEB1s);
title("TEB en fonction de (Eb/N0) (dB)");
xlabel("(Eb/N0) (dB)");
ylabel("TEB");
legend("TEB simulé chaine 1","TEB simulé chaine 3");

%% Ajout d'un filtre passe bas

% Création filtre passe bas
fc = 1500;
N = 101;
k = (-(N-1)/2 : (N-1)/2);
filtre_bas1 = 2*(fc/Fe)*sinc(2*k*(fc/Fe));

% Canal
r = filter(filtre_bas1, 1, x);

% Reception
z_bas = filter(hr, 1, r);

% Diagramme de l'oeil
eyediagram (z_bas(length(h):Nb*Ns), 2*Ns, 2*Ns);

%% Ajout d'un filtre passe haut

% Création filtre passe bas
fc = 3000;
filtre_bas2 = 2*(fc/Fe)*sinc(2*k*(fc/Fe));

% Canal
r = filter(filtre_bas2, 1, x);

% Reception
z_bas = filter(hr, 1, r);

% Diagramme de l'oeil
eyediagram (z_bas(length(h):Nb*Ns), 2*Ns, 2*Ns);






















