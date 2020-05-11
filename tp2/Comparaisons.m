%% Comparaison de modulatins sur fréquence proteuse

%% Initialisation
Eb_sur_N0_dB = [1:6];                       % vecteur Eb/N0 en dB
Eb_sur_N0 = 10.^(Eb_sur_N0_dB./10);         % vecteur Eb/N0 linéaire

Nerrlimite = 100;                           % nombre d'erreur limite pour une precision de 10%

Nbits = 1000;                               % Nombre de bits emis

Fe=144000;                                   % frequence d'echantillonnage
Rb=48000;                                   % debit binaire
retard_Ts=6;                                % Longueur du filtre en SRRCF en nombre de periode Ts
alphaSRRCF=0.5;                             % Coefficient roll off


%% Chaine 4-ASK

% Initialisation
TEB1s = zeros(1,length(Eb_sur_N0));         % Initialisation du TEB
Nerr=zeros(1,length(Eb_sur_N0_dB));         % Initialisation du nombre d'erreur
Sigma2 = zeros(1,length(Eb_sur_N0_dB));     % Initialisation de Sigma2
nbEssais=zeros(1,length(Sigma2));           % Initialisation du nombre d'essais
bits_estimes = zeros(1,Nbits);              % Initialisation de bits_estime
M = 4;                                      % Taille de la constellation
Rs = Rb/log2(M);                            % debit symbole
Ns=Fe/Rs;                                   % Nombre d'echantillons par periode porteuse
h=rcosdesign(alphaSRRCF,retard_Ts,Ns);      % Filtre d'emission
hr = fliplr(h);                             % Filtre de réception adapte
t0 = retard_Ts*Ns+1;                        % Instant de prise de décision pour d0 emis a t=0

% Emission
bits=randi([0 1],1,Nbits);      % Generation des bits
symboles = (2*bi2de(reshape(bits, 2, length(bits)/2).')-3).';
peigne= kron(symboles, [1, zeros(1,Ns-1)]); % Symboles du peigne de Dirac
peigne_allonge=[peigne zeros(1,retard_Ts*Ns)];
xe=filter(h,1,peigne_allonge); % Generation du signal passe-bas

% Canal
signal_canal = xe; % Signal du canal

% Reception
z=filter(hr,1,signal_canal);                               % Signal canal filtré
% Echantilonage
ze = z(t0:Ns:(Nbits/log2(M)+retard_Ts)*Ns);
        
% Decision
symbole_estimes = 3*(ze>2)-3*(ze<-2)+(ze>0 & ze<2)-(ze<0 & ze>-2);

% Demapping
bits_estimes = reshape(de2bi((symbole_estimes + 3)/2).',1,length(bits));        
        
TEB1 = sum(bits ~= bits_estimes); % Comparaison des bits
fprintf("Le TEB sans bruit vaut : %d \n", TEB1);




%% Chaine QPSK

% Initialisation
TEB2s = zeros(1,length(Eb_sur_N0));         % Initialisation du TEB
Nerr=zeros(1,length(Eb_sur_N0_dB));         % Initialisation du nombre d'erreur
Sigma2 = zeros(1,length(Eb_sur_N0_dB));     % Initialisation de Sigma2
nbEssais=zeros(1,length(Sigma2));           % Initialisation du nombre d'essais
bits_estimes = zeros(1,Nbits);              % Initialisation de bits_estime
M = 4;                                      % Taille de la constellation
Rs = Rb/log2(M);                            % debit symbole
Ns=Fe/Rs;                                   % Nombre d'echantillons par periode porteuse
h=rcosdesign(alphaSRRCF,retard_Ts,Ns);      % Filtre d'emission
hr = fliplr(h);                             % Filtre de réception adapte
t0 = retard_Ts*Ns+1;                        % Instant de prise de décision pour d0 emis a t=0

% Emission
bits=randi([0 1],1,Nbits);      % Generation des bits
symbI=2*bits(1:2:Nbits)-1;      % Codage symboles sur la voie I
symbQ=2*bits(2:2:Nbits)-1;      % Codage symboles sur la voie Q
symboles = symbI+1i*symbQ;
peigne= kron(symboles, [1, zeros(1,Ns-1)]); % Symboles du peigne de Dirac
peigne_allonge=[peigne zeros(1,retard_Ts*Ns)];
xe=filter(h,1,peigne_allonge); % Generation du signal passe-bas

% Canal
signal_canal = xe; % Signal du canal

% Reception
z=filter(hr,1,signal_canal);                               % Signal canal filtré
% Echantilonage
ze = z(t0:Ns:(Nbits/log2(M)+retard_Ts)*Ns);

% Demapping
% Les bits impaires correspondent à la partie réelle des symboles
bits_estimes(1:2:end) = real(ze) > 0;
% Les bits paires correspondent à la partie complexe des symboles
bits_estimes(2:2:end) = imag(ze) > 0;      
        
TEB2 = sum(bits ~= bits_estimes); % Comparaison des bits
fprintf("Le TEB sans bruit vaut : %d \n", TEB2);

%% chaine 8-PSK

% Initialisation
TEB3s = zeros(1,length(Eb_sur_N0));         % Initialisation du TEB
Nerr=zeros(1,length(Eb_sur_N0_dB));         % Initialisation du nombre d'erreur
Sigma2 = zeros(1,length(Eb_sur_N0_dB));     % Initialisation de Sigma2
nbEssais=zeros(1,length(Sigma2));           % Initialisation du nombre d'essais
bits_estimes = zeros(1,Nbits);              % Initialisation de bits_estime
M = 8;                                      % Taille de la constellation
Rs = Rb/log2(M);                            % debit symbole
Ns=Fe/Rs;                                   % Nombre d'echantillons par periode porteuse
h=rcosdesign(alphaSRRCF,retard_Ts,Ns);      % Filtre d'emission
hr = fliplr(h);                             % Filtre de réception adapte
t0 = retard_Ts*Ns+1;                        % Instant de prise de décision pour d0 emis a t=0

% Emission
data=randi([0 M-1],1,Nbits);      % Generation des bits
symboles = pskmod(bits,M,pi/M);
peigne= kron(symboles, [1, zeros(1,Ns-1)]); % Symboles du peigne de Dirac
peigne_allonge=[peigne zeros(1,retard_Ts*Ns)];
xe=filter(h,1,peigne_allonge); % Generation du signal passe-bas

% Canal
signal_canal = xe; % Signal du canal

% Reception
z=filter(hr,1,signal_canal);                               % Signal canal filtré
% Echantilonage
ze = z(t0:Ns:(Nbits/log2(M)+retard_Ts)*Ns);

% Demapping
symboles_estimes = pskdemod(ze,M,pi/M);        
        
TEB3 = sum(symboles ~= symboles_estimes); % Comparaison des bits
fprintf("Le TEB sans bruit vaut : %d \n", TEB3);

%% chaine 16-QAM

% Initialisation
TEB4s = zeros(1,length(Eb_sur_N0));         % Initialisation du TEB
Nerr=zeros(1,length(Eb_sur_N0_dB));         % Initialisation du nombre d'erreur
Sigma2 = zeros(1,length(Eb_sur_N0_dB));     % Initialisation de Sigma2
nbEssais=zeros(1,length(Sigma2));           % Initialisation du nombre d'essais
bits_estimes = zeros(1,Nbits);              % Initialisation de bits_estime
M = 16;                                     % Taille de la constellation
Rs = Rb/log2(M);                            % debit symbole
Ns=Fe/Rs;                                   % Nombre d'echantillons par periode porteuse
h=rcosdesign(alphaSRRCF,retard_Ts,Ns);      % Filtre d'emission
hr = fliplr(h);                             % Filtre de réception adapte
t0 = retard_Ts*Ns+1;                        % Instant de prise de décision pour d0 emis a t=0

% Emission
bits=randi([0 M-1],1,Nbits);      % Generation des bits
symboles = qammod(bits,M,pi/M);
peigne= kron(symboles, [1, zeros(1,Ns-1)]); % Symboles du peigne de Dirac
peigne_allonge=[peigne zeros(1,retard_Ts*Ns)];
xe=filter(h,1,peigne_allonge); % Generation du signal passe-bas

% Canal
signal_canal = xe; % Signal du canal

% Reception
z=filter(hr,1,signal_canal);                               % Signal canal filtré
% Echantilonage
ze = z(t0:Ns:(Nbits/log2(M)+retard_Ts)*Ns);

% Demapping
symboles_estimes = qamdemod(ze,M,pi/M);        
        
TEB4 = sum(symboles ~= symboles_estimes); % Comparaison des bits
fprintf("Le TEB sans bruit vaut : %d \n", TEB4);
