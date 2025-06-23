% Confronto tra dinamica esatta e integrazione numerica 
% dell'equazione di Landau-Lifshitz-Gilbert (LLG)
% Obiettivo: valutare durante la dinamica l'errore commesso su Sx e Sz
%            a seguito dell'integrazione numerica usando il metodo di Heun          
% -----------------------------------------------------------------------

% Definizione della funzione per implementare la matrice M(S)
% tramite interpretazione del prodotto vettoriale in termini di 
% matrice 3x3
function M_S = M_S(S, gamma, alpha)

    % S_cross: matrice anti‑simmetrica che realizza [S]_x
    S_cross = [   0    -S(3)   S(2);
             S(3)      0    -S(1);
            -S(2)    S(1)     0  ];

    % costante prefattore
    kappa = gamma / (1 + alpha^2);

    % costruzione di M(S)
    M_S = - kappa * ( S_cross + alpha * (S_cross * S_cross) );
    
end

% Parametri fisici
gamma = 1.76e11;    % fattore giromagnetico [rad/(s T)]
alpha = 0.1;        % parametro di smorzamento [adimensionato]
H = [0; 0; 10];     % campo magnetico applicato in [T] lungo z
mu = 10e-24;        % magnetone di Bohr [J/T]

% Parametri di integrazione
dt = 1e-15;           % durata simulazione [s]
t_max = 50e-12;       % passo temporale [s]
t = 0:dt:t_max;       % vettore dei tempi
N = length(t) - 1;    % numero di passi temporali

% Preallocazioni
S = zeros(3, N);          % traiettoria numerica 
S_exact = zeros(3, N);    % soluzione esatta
err_x = zeros(1, N);      % errore su Sx
err_z = zeros(1, N);      % errore su Sz

% Condizione iniziale
S(:,1) = [1; 0; 0];
S_exact(:,1)  = [1; 0; 0];

% Funzione sech
sech = @(x) 1 ./ cosh(x);

% Costanti per soluzione esatta
denom = 1 + alpha^2;
lambda = (alpha * gamma * H(3)) / denom;
omega = (gamma * H(3)) / denom;

% Integrazione numerica (Heun stocastico con σ=0)
for n = 1:N
    
    % Costruzione del sistema al passo n-esimo
    M_n = M_S( S(:,n), gamma, alpha );  % M(S_n)
    drift  = M_n * H;                        

    % Step predittore
    S_tilde = S(:,n) + drift * dt;

    % Step correttore
    M_tilde = M_S( S_tilde, gamma, alpha ); % M(S_tilde)
    drift_tilde = M_tilde * H;

    % Step complessivo secondo Heun stocastico
    S(:,n+1) = S(:,n) + 0.5*(drift + drift_tilde)*dt;

    % Normalizzazione
    S(:,n+1) = S(:,n+1) / norm(S(:,n+1));

    % Soluzione esatta al tempo t(n+1)
    tau = t(n+1);
    S_exact(1,n+1) = sech(lambda * tau) * cos(omega * tau);
    S_exact(2,n+1) = sech(lambda * tau) * sin(omega * tau);
    S_exact(3,n+1) = tanh(lambda * tau);

    % Errori componenti
    err_x(n+1) = S(1,n+1) - S_exact(1,n+1);
    err_z(n+1) = S(3,n+1) - S_exact(3,n+1);
end

% Grafico errori sulle componenti di spin
figure;
plot(t * 1e12, [err_x; err_z]', 'LineWidth', 2);
ax = gca; set(ax, 'FontSize', 14);
xlabel('Tempo [ps]', 'FontSize', 16);
ylabel('Errore componenti di spin', 'FontSize', 16);
legend({'\Delta S_x','\Delta S_z'}, 'FontSize', 14, 'Location','best');
grid on; ax.GridAlpha = 0.3;