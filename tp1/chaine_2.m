% Deuxième chaine étudiée : impact du choix de filtre de réception
close all
clear all

%% Initialisation des constantes
Nb = 10000; % Nombre de bits
Ns = 8;    % Nombre d'échantillon par période symbole
h = ones(1, Ns); % Répertoire impulsionnelle du filtre de mise en forme
hr1 = fliplr(h);  %Filtre de réception adapté chaine 1

hr2 = [ones(1, Ns/2), zeros(1, Ns/2)];  %Filtre de réception adapté chaine 2
t01 = length(h);    % Cf chaine 1
t02 = linspace(Ns/2, Ns, Nb);   % t0=[Ns/2 Ns]
                                % t0 represente l'instant de prise de décision
                                % pour le symbole a0 emis a t=0

hr2 = [ones(1, Ns/2), zeros(1, Ns/2)];  %Filtre de réception adapté chaine 2
t01 = length(h);
t02 = linspace(Ns/2, Ns, Nb);

%% Génération des bits et Mapping
bits = randi([0,1],1,Nb);
symboles = 2*bits - 1;

peigne_dirac = kron(symboles, [1, zeros(1,Ns-1)]);
x = filter(h, 1, peigne_dirac);


%% Densité spectrale de puissance
X2 = fft(x);
DSP2 = 1/(Nb*Ns) * abs(X2).^2;

X1 = fft(x);
DSP1 = 1/(Nb*Ns) * abs(X1).^2;

figure; 
semilogy(linspace(-0.5, 0.5, length(DSP1)),fftshift(DSP1)); hold on;
semilogy(linspace(-0.5, 0.5, length(DSP2)),fftshift(DSP2));
title("DSP du signal transmis");
xlabel("Fréquence normalisée");
ylabel("DSP(f)");
legend("DSP chaine de référence","DSP chaine étudiée");

%% Canal
r = x;

%% Reception
z = filter(hr2, 1, r);
z1 = filter(hr1, 1, r);


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
Eb_sur_N0_dB = linspace(0,6,50);            % Mise en place du rapport signal
Eb_sur_N0 = 10.^(Eb_sur_N0_dB./10);         % sur bruit variat de 1 à 6
TEB1s = zeros(1,length(Eb_sur_N0));         % Initialisation du TEB de la chaine 1
TEB2s = zeros(1,length(Eb_sur_N0));         % Initialisation du TEB de la chaine 2

M=2;                                        % Nombre de symboles
Pr = mean(abs(x).^2);                       % Calcul de la puissance du signal
Sigma2 = Pr*Ns./(2*log2(M)*Eb_sur_N0);      % Calcul de la variance du bruit
% si X est N(0,1)
% alors aX+b est N(b,a^2)


nbEssais1=zeros(1,length(Sigma2)); % Initialisation du nombre d'essais
Nerr1=zeros(1,length(Sigma2));     % Initialisation du nombre d'erreur


nbEssais2=zeros(1,length(Sigma2)); % Initialisation du nombre d'essais
Nerr2=zeros(1,length(Sigma2));     % Initialisation du nombre d'erreur

Nelimite = 1000;
% TEB de la chaîne de référence
for i = 1:length(Sigma2)
    
    while (Nerr1(i) < Nelimite)
        % Canal avec bruit AWGN
        r = x + sqrt(Sigma2(i))*randn(1,length(x)); 
        % Réception
        z = filter(hr1, 1, r); 
        % Echantilonage
        ze = z(t01:Ns:Ns*Nb); 
        % Décision
        bits_estimes = (ze > 0);
        NerrActuel = sum(bits ~= bits_estimes);
        Nerr1(i) = Nerr1(i) + NerrActuel;
        nbEssais1(i) = nbEssais1(i) + 1;
    end
    TEB1s(i) = Nerr1(i)/(nbEssais1(i)*Nb);
end
% TEB de la chaîne 2
for i = 1:length(Sigma2)

    while (Nerr2(i) < Nelimite)
        % Canal avec bruit AWGN
        r = x + sqrt(Sigma2(i))*randn(1,length(x)); 
        % Réception
        z = filter(hr2, 1, r); 
        % Echantilonage
        ze = z(t02:Ns:Ns*Nb); 
        % Décision
        bits_estimes = (ze > 0);
        NerrActuel = sum(bits ~= bits_estimes);
        Nerr2(i) = Nerr2(i) + NerrActuel;
        nbEssais2(i) = nbEssais2(i) + 1;
    end
    TEB2s(i) = Nerr2(i)/(nbEssais2(i)*Nb);
end

TEB_theo = qfunc(sqrt(Eb_sur_N0));
Pb = TEB_theo; % Puissance du bruit
variance_simu=Pb.*(1-Pb)./(nbEssais2.*Nb); % Variance de la simu

figure;
semilogy(Eb_sur_N0_dB,TEB2s,'r+'); hold on;
semilogy(Eb_sur_N0_dB,TEB_theo,'g');
semilogy(Eb_sur_N0_dB,Pb+sqrt(variance_simu),'c')
semilogy(Eb_sur_N0_dB,Pb-sqrt(variance_simu),'c')
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
