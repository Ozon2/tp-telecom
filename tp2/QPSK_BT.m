% QPSK en bande transposee

clear all
close all
clc

% Initialisation
Eb_sur_N0_dB = [1:6]; % vecteur Eb/N0 en dB
Eb_sur_N0 = 10.^(Eb_sur_N0_dB./10);

Nerr=zeros(1,length(Eb_sur_N0_dB)); % Initialisation du nombre d'erreur
Nerrlimite = 100; % nombre d'erreur limite pour une precision de 10%

Nbits = 10 ; % Nombre de bits emis
M = 4 ; % Taille de la constellation
offsetQPSK=pi/4;

Fe=10000; % frequence d'echantillonnage
Rs=1000; % debit symbole
Ns=Fe/Rs; % Nombre d'echantillons par periode porteuse
retard_Ts=6; % Longueur du filtre en SRRCF en nombre de periode Ts
alphaSRRCF=0.35; % Coefficient roll off
h=rcosdesign(alphaSRRCF,retard_Ts,Ns); % Filtre d'emission
hr = fliplr(h); % Filtre de reception adapte
t0 = retard_Ts*Ns+1; % t0 represente l'instant de prise de décision
                     % pour le symbole a0 emis a t=0

fp=2000; % frequence porteuse
fpN=fp/Fe;    % frequence porteuse normalisee
Sigma2 = zeros(1,length(Eb_sur_N0_dB)); % Initialisation de Sigma2

for k=1:length(Eb_sur_N0_dB)
    
    %while (Nerr(k)<Nerrlimite)

        % Emission
        bits=randi([0 1],1,Nbits);  % Generation des bits
        symbI=1-2*bits(1:2:Nbits); % Codage symboles sur la voie I
        symbQ=1-2*bits(2:2:Nbits); % Codage symboles sur la voie Q
        symbole = symbI+1i*symbQ;
        peigne=zeros(1,Nbits/log2(M)*Ns);
        peigne(1:Ns:Nbits/log2(M)*Ns)=symbole; % Symboles du peigne de Dirac
        signal_passe_bas=filter(h,1,peigne); % Generation du signal passe-bas
        
        %Transposition sur fréquence porteuse
        signal_passe_bande=real(exp(1i*2*pi*fpN*[1:Nbits/log2(M)*Ns]).*signal_passe_bas);

        % Canal
        Pr = mean(abs(signal_passe_bande).^2);      % Calcul de la puissance du signal
        Sigma2(k) = Pr*Ns./(2*log2(M)*Eb_sur_N0(k));      % Calcul de la variance du bruit
        signal_canal = signal_passe_bande + sqrt(Sigma2(k))*randn(1,length(signal_passe_bas)); % Signal du canal
        
        % Retour en Bande de base
        signal_reel = signal_canal .* cos(2*pi*fp*[1:Nbits/log2(M)*Ns]);
        signal_complexe = signal_canal .* sin(2*pi*fp*[1:Nbits/log2(M)*Ns]);
        
        % Reception
        z = filter(hr,1,signal_reel)-1i*filter(hr,1,signal_reel); % Signal recu
        
        % Echantilonage
        ze = z(t0:Ns:Ns*Nbits/2);
        
        % Decision
        symbole_estime = (real(ze) > 0) + 1i*(imag(ze) > 0) - (real(ze) < 0) - 1i*(imag(ze) < 0);
        bits_estime = 
    %end
end


nfft=2^nextpow2(Nbits/log2(M)*Ns); % Nombre de points de la FFT
w = window(@blackmanharris, length(signal_passe_bande)); % Fenetre du periodogramme
[Px1,F1] = periodogram(signal_passe_bas,w,nfft); % periodogramme du signal passe bas
[Px2,F2] = periodogram(signal_passe_bande,w,nfft); % periodogramme du signal passe bande
semilogy(F1, fftshift(Px1));
hold on
semilogy(F2, fftshift(Px2), 'r');
xlabel('frequence normalisee')
grid
legend('PSD signal passe-bas','PSD signal passe-bande')
