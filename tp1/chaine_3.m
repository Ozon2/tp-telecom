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

% Densité spectrale de puissance chaine 1
X = fft(x);
DSP1 = 1/(Nb*Ns) * abs(X).^2;

% Canal chaine 1
r = x;

% Reception chaine 1
z = filter(hr, 1, r);


% TEB avec bruit chaine 1
Eb_sur_N0_dB = linspace(0,6,50);
Eb_sur_N0 = 10.^(Eb_sur_N0_dB./10);
TEB1s = zeros(1,length(Eb_sur_N0));
Pr = mean(abs(x).^2);
Sigma2 = Pr*Ns./(2*Eb_sur_N0);

Nelimite = 1000;
for i = 1:length(Sigma2)
    Nerr = 0;
    nbEssais = 0;
    while (Nerr < Nelimite)
        % Canal avec bruit AWGN
        r = x + sqrt(Sigma2(i))*randn(1,length(x));
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
t0 = SPAN*Ns+1; % t0 represente l'instant de prise de décision
                % pour le symbole a0 emis a t=0

%% Génération des bits et Mapping
bits = randi([0,1],1,Nb);
symboles = 2*bits - 1;

peigne_dirac = kron(symboles, [1, zeros(1,Ns-1)]);
x = filter(h, 1, peigne_dirac);

%% Densité spectrale de puissance
X = fft(x);
DSP = 1/(Nb*Ns) * abs(X).^2; 
figure; semilogy(linspace(-0.5, 0.5, length(DSP)),fftshift(DSP)); hold on;
semilogy(linspace(-0.5, 0.5, length(DSP1)),fftshift(DSP1));
title("DSP du signal transmis");
xlabel("Fréquence normalisée");
ylabel("DSP(f)");
legend("DSP chaine étudiée","DSP chaine de référence");

%% Canal
r = x;

%% Reception
z = filter(hr, 1, r);

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
Eb_sur_N0_dB = linspace(0,6,50);            % Mise en place du rapport signal
Eb_sur_N0 = 10.^(Eb_sur_N0_dB./10);         % sur bruit variat de 1 à 6
TEBs = zeros(1,length(Eb_sur_N0));          % Initialisation du TEB

M=2;                                        % Nombre de symboles
Pr = mean(abs(x).^2);                       % Calcul de la puissance du signal
Sigma2 = Pr*Ns./(2*log2(M)*Eb_sur_N0);      % Calcul de la variance du bruit
% si X est N(0,1)
% alors aX+b est N(b,a^2)

nbEssais=zeros(1,length(Sigma2)); % Initialisation du nombre d'essais
Nerr=zeros(1,length(Sigma2));     % Initialisation du nombre d'erreur

Nelimite = 100;  % Nombre d'erreurs attendues pour une précision de 
                % e=2/sqrt(Nelimite) 
for i = 1:length(Sigma2)

    while (Nerr(i) < Nelimite)
        % Canal avec bruit AWGN
        r = x + sqrt(Sigma2(i))*randn(1,length(x)); 
        % Réception
        z = filter(hr, 1, r); 
        % Echantilonage
        ze = z(t0:Ns:Ns*Nb); 
        % Décision
        bits_estimes = (ze > 0);
        NerrActuel = sum(bits(1:Nb-SPAN) ~= bits_estimes);
        Nerr(i) = Nerr(i) + NerrActuel;
        nbEssais(i) = nbEssais(i) + 1;
    end
    TEBs(i) = Nerr(i)/(nbEssais(i)*Nb);
end

TEB_theo = qfunc(sqrt(2*Eb_sur_N0));
Pb = TEB_theo; % Puissance du bruit
variance_simu=Pb.*(1-Pb)./(nbEssais.*Nb); % Variance de la simu

figure;
semilogy(Eb_sur_N0_dB,TEBs,'r+'); hold on;
semilogy(Eb_sur_N0_dB,TEB_theo,'g');
semilogy(Eb_sur_N0_dB,Pb+sqrt(variance_simu),'c')
semilogy(Eb_sur_N0_dB,Pb-sqrt(variance_simu),'c')
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

%% Ajout d'un filtre passe bas de fréquence de coupure 1500Hz

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
title("Diagramme de l'oeil avec un filtre passe bas de fréquence de coupure 1500 Hz");

%% Ajout d'un filtre passe bas de fréquence de coupure 3000Hz

% Création filtre passe bas
fc = 3000;
filtre_bas2 = 2*(fc/Fe)*sinc(2*k*(fc/Fe));

% Canal
r = filter(filtre_bas2, 1, x);

% Reception
z_bas = filter(hr, 1, r);

% Diagramme de l'oeil
eyediagram (z_bas(length(h):Nb*Ns), 2*Ns, 2*Ns);
title("Diagramme de l'oeil avec un filtre passe bas de fréquence de coupure 3000 Hz");





















