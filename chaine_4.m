% Quatrième chaine étudiée : impacte du choix du mapping
close all
clear all

%% Initialisation des constantes
Nb = 10000; % Nombre de bits
Ns = 8;    % Nombre d'échantillon par période symbole
h = ones(1, Ns); % Répertoire impulsionNslle du filtre de mise en forme
hr = fliplr(h);  %Filtre de réception adapté
t0 = length(h);

%% Première chaine

% Génération des bits et Mapping chaine 1
bits = randi([0,1],1,Nb);
symboles = 3*(2*bits - 1);
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
Eb_sur_N0_dB = linspace(0,6,100);
Eb_sur_N0 = 10.^(Eb_sur_N0_dB./10);
TEB1s = zeros(1,length(Eb_sur_N0));
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


%% Génération des bits et Mapping
bits = randi([0,1],1,Nb);
symboles = (2*bi2de(reshape(bits, 2, length(bits)/2).')-3).';

peigNs_dirac = kron(symboles, [1, zeros(1,Ns-1)]);
x = filter(h, 1, peigNs_dirac);

%% Canal
r = x;

%% Reception
z = filter(hr, 1, r);


%% Densité spectrale de puissance
Z = fft(z);
DSP = 1/(Nb*Ns) * abs(Z).^2; 
figure;
semilogy(linspace(-0.5, 0.5, length(DSP1)),fftshift(DSP1)); hold on;
semilogy(linspace(-0.5, 0.5, length(DSP)),fftshift(DSP));
title("DSP du signal transmis");
xlabel("Fréquence normalisée");
ylabel("DSP(f)");
legend("DSP chaine étudiée","DSP chaine de référence");

%% Echantillonnage
ze = z(t0:Ns:Ns*Nb/2);

%% Tracé de la chaine
figure;
plot((1:Nb*Ns/2)/Ns, peigNs_dirac,'LineWidth',2);
hold on
plot((1:Nb*Ns/2)/Ns, x, 'r--','LineWidth',2);
plot((1:Nb*Ns/2)/Ns, z, 'b.-','LineWidth',2);
plot((t0:Ns:Nb*Ns/2)/Ns, ze ,'gp','LineWidth',3);
legend ("PeigNs de Dirac", "x(t)", "z(t)", "ze(t)");
xlabel("Temps (en période T_s)");
title("chaine de référence");
grid

%% Diagramme de l'oeil
eyediagram (z(length(h):Nb*Ns/2)/8, 2*Ns, 2*Ns);

%% Décision
symbole_estimes = 3*(ze>16)-3*(ze<-16)+(ze>0 & ze<16)-(ze<0 & ze>-16);
bits_estimes = reshape(de2bi((symbole_estimes + 3)/2).',1,length(bits));
TEB = sum(bits ~= bits_estimes)/Nb;
fprintf("Le TEB sans bruit vaut : %d \n", TEB);

%% TEB avec bruit
Eb_sur_N0_dB = linspace(0,6,100);
Eb_sur_N0 = 10.^(Eb_sur_N0_dB./10);
TESs = zeros(1,length(Eb_sur_N0));
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
         z = filter(hr, 1, r);
        % Echantilonage
        ze = z(t0:Ns:Ns*Nb/2);
        % Décision
        symbole_estimes = 3*(ze>16)-3*(ze<-16)+(ze>0 & ze<16)-(ze<0 & ze>-16);
        bits_estimes = reshape(de2bi((symbole_estimes + 3)/2).',1,length(bits));
        NerrActuel = sum(bits ~= bits_estimes);
        Nerr = Nerr + NerrActuel;
        nbEssais = nbEssais + 1;
    end
    TESs(i) = Nerr/(nbEssais*Nb);
end
TEBs = TESs/log2(4);

TES_theo = 2*(3/4)*qfunc(sqrt(6*log2(4)/(4^2-1)*Eb_sur_N0));
TEB_theo = TES_theo/log2(4);

figure;
semilogy(Eb_sur_N0_dB,TESs); hold on;
semilogy(Eb_sur_N0_dB,TES_theo,'r*');
title("TES en fonction de (Eb/N0) (dB)");
legend("TES expérimental","TES théorique");
xlabel("(Eb/N0) (dB)");
ylabel("TES");
legend("TES simulé","TES théorique");

figure;
semilogy(Eb_sur_N0_dB,TEBs); hold on;
semilogy(Eb_sur_N0_dB,TEB_theo,'r*');
title("TEB en fonction de (Eb/N0) (dB)");
legend("TEB expérimental","TES théorique");
xlabel("(Eb/N0) (dB)");
ylabel("TEB");
legend("TEB simulé","TEB théorique");

figure;
semilogy(Eb_sur_N0_dB,TEBs); hold on;
semilogy(Eb_sur_N0_dB,TEB1s);
title("TEB en fonction de (Eb/N0) (dB)");
xlabel("(Eb/N0) (dB)");
ylabel("TEB");
legend("TEB simulé chaine 1","TEB simulé chaine 4");
