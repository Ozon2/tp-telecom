%% QPSK en bande transposee

clear all
close all
clc

% Initialisation
Eb_sur_N0_dB = [1:6];                       % vecteur Eb/N0 en dB
Eb_sur_N0 = 10.^(Eb_sur_N0_dB./10);         % vecteur Eb/N0 linéaire
TEB1s = zeros(1,length(Eb_sur_N0));          % Initialisation du TEB


Nerr=zeros(1,length(Eb_sur_N0_dB));         % Initialisation du nombre d'erreur
Nerrlimite = 100;                           % nombre d'erreur limite pour une precision de 10%

Nbits = 1000;                                % Nombre de bits emis
M = 4 ;                                     % Taille de la constellation

Fe=10000;                                   % frequence d'echantillonnage
Rs=1000;                                    % debit symbole
Ns=Fe/Rs;                                   % Nombre d'echantillons par periode porteuse
retard_Ts=6;                                % Longueur du filtre en SRRCF en nombre de periode Ts
alphaSRRCF=0.35;                            % Coefficient roll off
h=rcosdesign(alphaSRRCF,retard_Ts,Ns);      % Filtre d'emission
hr = fliplr(h);                             % Filtre de réception adapte
t0 = retard_Ts*Ns+1;                        % Instant de prise de décision pour d0 emis a t=0
fp=2000;                                    % frequence porteuse
fpN=fp/Fe;                                  % frequence porteuse normalisee
Sigma2 = zeros(1,length(Eb_sur_N0_dB));     % Initialisation de Sigma2
nbEssais=zeros(1,length(Sigma2));           % Initialisation du nombre d'essais
bits_estimes = zeros(1,Nbits);              % Initialisation de bits_estime

temps = [1:(Nbits/log2(M)+retard_Ts)*Ns];

% Création filtre passe bas
fc = fp;
N = 101;
k = (-(N-1)/2 : (N-1)/2);
filtre_bas = 2*(fc/Fe)*sinc(2*k*(fc/Fe));


% Emission
bits=randi([0 1],1,Nbits);      % Generation des bits
symbI=2*bits(1:2:Nbits)-1;      % Codage symboles sur la voie I
symbQ=2*bits(2:2:Nbits)-1;      % Codage symboles sur la voie Q
symboles = symbI+1i*symbQ;
peigne= kron(symboles, [1, zeros(1,Ns-1)]); % Symboles du peigne de Dirac
peigne_allonge=[peigne zeros(1,retard_Ts*Ns)];
xe=filter(h,1,peigne_allonge); % Generation du signal passe-bas
        
%Transposition sur fréquence porteuse
x=real(exp(1i*2*pi*fpN*temps).*xe); % Generation du signal passe-bande

Pr = mean(abs(x).^2);                       % Calcul de la puissance du signal recu
Sigma2 = Pr*Ns./(2*log2(M)*Eb_sur_N0);      % Calcul de la variance du bruit

for k=1:length(Eb_sur_N0_dB)
    
    while (Nerr(k) < Nerrlimite)
        
        

        % Canal
        signal_canal = x + sqrt(Sigma2(k))*randn(1,length(x)); % Signal du canal
        
        % Retour en Bande de base
        signal_reel = signal_canal .* cos(2*pi*fpN*temps);     % correspond à x(t)*cos(2*pi*fp*t)
        signal_complexe = signal_canal .* sin(2*pi*fpN*temps); % correspond à x(t)*sin(2*pi*fp*t)
        
        signal_recu_reel = conv(signal_reel,filtre_bas,'same'); % Filtrage passe bas de la partie réelle
        signal_recu_complexe = conv(signal_complexe,filtre_bas,'same'); % Filtrage passe bas de la partie complexe
        
        % Reception
        signal_recu = signal_recu_reel - 1i*signal_recu_complexe; % signal recu
        z=filter(hr,1,signal_recu);                               % Signal recu filtré
        % Echantilonage
        ze = z(t0:Ns:(Nbits/log2(M)+retard_Ts)*Ns);
        
        % Decision
        symboles_estimes = (real(ze) > 0) - (real(ze) < 0) + 1i*((imag(ze) > 0) - (imag(ze) < 0));
        
        % Demapping
        % Les bits impaires correspondent à la partie réelle des symboles
        bits_estimes(1:2:end) = real(symboles_estimes) > 0;
        % Les bits paires correspondent à la partie complexe des symboles
        bits_estimes(2:2:end) = imag(symboles_estimes) > 0;
        
        
        NerrActuel = sum(bits ~= bits_estimes); % Comparaison des bits
        Nerr(k) = Nerr(k) + NerrActuel; % Calcul du nombre d'erreur
        
        nbEssais(k) = nbEssais(k) + 1; % Calcul du nombre d'essais
    end
    TEB1s(k) = Nerr(k)/(nbEssais(k)*Nbits);
end


nfft=2^nextpow2(Nbits/log2(M)*Ns);      % Nombre de points de la FFT
axe_f=linspace(-0.5,0.5,nfft) ;         % Axe des frequences
w = window(@blackmanharris, length(x)); % Fenetre du periodogramme
[Px1,F1] = periodogram(xe,w,nfft);      % periodogramme du signal passe bas
[Px2,F2] = periodogram(x + eps*1i,w,nfft);       % periodogramme du signal passe bande

figure
plot(ze,'r+'); hold on
plot([-1.5,1.5], [0,0], 'g-');
plot([0,0],[-1.5,1.5],'g-');
title("Constellation");
xlabel("Partie réelle");
ylabel("Partie imaginaire");
legend("symboles");
axis([-1.5 1.5 -1.5 1.5]);

figure
x_r = real(xe);
subplot(2,1,1); plot(x_r);
title("Signal en phase");
x_c = imag(xe);
subplot(2,1,2); plot(x_c);
title("Signal en quadrature de phase");

figure
plot(x);
title("Signal sur fréquence porteuse");

figure
semilogy(axe_f,fftshift(abs(fft(xe,nfft)).^2))
hold on
semilogy(axe_f,fftshift(abs(fft(x,nfft)).^2),'r')
xlabel('frequence normalisee')
grid
legend('DSP signal passe-bas','DSP signal passe-bande')

figure
semilogy(F1, fftshift(Px1));
hold on
semilogy(F2, fftshift(Px2), 'r');
title("Periodogramme de signaux");
xlabel('frequence normalisee');
grid on
legend('Periodogramme signal passe-bas','Periodogramme signal passe-bande');


TEB_theo = qfunc(sqrt(2*Eb_sur_N0));
Pb = TEB_theo; % Puissance du bruit
variance_simu=Pb.*(1-Pb)./(nbEssais.*Nbits); % Variance de la simu

figure;
semilogy(Eb_sur_N0_dB,TEB1s,'r+'); hold on;
semilogy(Eb_sur_N0_dB,TEB_theo,'g');
semilogy(Eb_sur_N0_dB,Pb+sqrt(variance_simu),'c')
semilogy(Eb_sur_N0_dB,Pb-sqrt(variance_simu),'c')
title("TEB en fonction de (Eb/N0) (dB)");
xlabel("(Eb/N0) (dB)");
ylabel("TEB");
legend("TEB simulé","TEB théorique");


%% QPSK en passe-bas équivalent

% Initialisation
TEB2s = zeros(1,length(Eb_sur_N0));          % Initialisation du TEB

% Emission
% On réutilise les memes bits que pour la chaine 1
symbI=2*bits(1:2:Nbits)-1;      % Codage symboles sur la voie I
symbQ=2*bits(2:2:Nbits)-1;      % Codage symboles sur la voie Q
symboles = symbI+1i*symbQ;
peigne= kron(symboles, [1, zeros(1,Ns-1)]); % Symboles du peigne de Dirac
peigne_allonge=[peigne zeros(1,retard_Ts*Ns)];
xe=filter(h,1,peigne_allonge); % Generation du signal passe-bas

Pr = mean(abs(x).^2);                       % Calcul de la puissance du signal recu
Sigma2 = Pr*Ns./(2*log2(M)*Eb_sur_N0);      % Calcul de la variance du bruit

for k=1:length(Eb_sur_N0_dB)
    
    while (Nerr(k) < Nerrlimite)
        
        

        % Canal
        signal_canal = xe + sqrt(Sigma2(k))*(randn(1,length(xe)) + 1i*randn(1,length(xe))); % Signal du canal

        % Reception
        z=filter(hr,1,signal_canal);                               % Signal canal filtré
        % Echantilonage
        ze = z(t0:Ns:(Nbits/log2(M)+retard_Ts)*Ns);
        
        % Decision
        symboles_estimes = (real(ze) > 0) - (real(ze) < 0) + 1i*((imag(ze) > 0) - (imag(ze) < 0));
        
        % Demapping
        % Les bits impaires correspondent à la partie réelle des symboles
        bits_estimes(1:2:end) = real(symboles_estimes) > 0;
        % Les bits paires correspondent à la partie complexe des symboles
        bits_estimes(2:2:end) = imag(symboles_estimes) > 0;
        
        
        NerrActuel = sum(bits ~= bits_estimes); % Comparaison des bits
        Nerr(k) = Nerr(k) + NerrActuel; % Calcul du nombre d'erreur
        
        nbEssais(k) = nbEssais(k) + 1; % Calcul du nombre d'essais
    end
    TEB2s(k) = Nerr(k)/(nbEssais(k)*Nbits);
end


nfft=2^nextpow2(Nbits/log2(M)*Ns);      % Nombre de points de la FFT
axe_f=linspace(-0.5,0.5,nfft) ;         % Axe des frequences
w = window(@blackmanharris, length(xe)); % Fenetre du periodogramme
[Px1,F1] = periodogram(imag(xe),w,nfft);      % periodogramme du signal passe bas

figure
plot(ze,'r+'); hold on
plot([-2,2], [0,0], 'g-');
plot([0,0],[-2,2],'g-');
title("Constellation");
xlabel("Partie réelle");
ylabel("Partie imaginaire");
legend("symboles");
axis([-2 2 -2 2]);

figure
x_r = real(xe);
subplot(2,1,1); plot(x_r);
title("Signal en phase");
x_c = imag(xe);
subplot(2,1,2); plot(x_c);
title("Signal en quadrature de phase");

figure
semilogy(axe_f,fftshift(abs(fft(xe,nfft)).^2))
hold on
semilogy(axe_f,fftshift(abs(fft(x,nfft)).^2),'r')
xlabel('frequence normalisee')
grid
legend('DSP enveloppe complexe passe-bas','DSP signal sur fréquence porteuse')

figure
semilogy(F1, fftshift(Px1));
hold on
semilogy(F2, fftshift(Px2), 'r');
title("Periodogramme de signaux");
xlabel('frequence normalisee');
grid on
legend('Periodogramme enveloppe complexe passe-bas','Periodogramme signal sur fréquence porteuse');


TEB_theo = qfunc(sqrt(2*Eb_sur_N0));
Pb = TEB_theo; % Puissance du bruit
variance_simu=Pb.*(1-Pb)./(nbEssais.*Nbits); % Variance de la simu

figure;
semilogy(Eb_sur_N0_dB,TEB1s,'b+'); hold on
semilogy(Eb_sur_N0_dB,TEB2s,'r+');
semilogy(Eb_sur_N0_dB,TEB_theo,'g');
semilogy(Eb_sur_N0_dB,Pb+sqrt(variance_simu),'c')
semilogy(Eb_sur_N0_dB,Pb-sqrt(variance_simu),'c')
title("TEB en fonction de (Eb/N0) (dB)");
xlabel("(Eb/N0) (dB)");
ylabel("TEB");
legend("TEB simulé fréquence porteuse","TEB simulé passe-bas equivalent","TEB théorique");
