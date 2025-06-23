%--------------------------------------------------------------------------
% Simulazione dell'equazione stocastica di Landau-Lifshitz-Gilbert (LLG)
% Metodo: Heun stocastico (predittore-correttore)
% Obiettivo: simulare la dinamica di uno spin classico in presenza di 
%            campo magnetico e fluttuazioni termiche (rumore bianco)
%--------------------------------------------------------------------------

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
T = 0.1;            % temperatura [K]
kB = 1.38e-23;      % costante di Boltzmann [J/K]
mu = 9.274e-24;     % magnetone di Bohr [J/T]

% intensità rumore termico (rumore bianco)
sigma = sqrt( 2 * alpha * kB * T / (gamma * mu * (1 + alpha^2)) ) 

% Parametri di integrazione
dt = 1e-15;           % durata simulazione [s]
t_max = 50e-12;       % passo temporale [s]
t = 0:dt:t_max;       % vettore dei tempi
N = length(t) - 1;    % numero di step

% Preallocazione
S = zeros(3, N);    % traiettoria di spin

% Condizione iniziale
S(:,1) = [1; 0; 0]; % spin inizialmente lungo x

% costante 
sqrt_dt = sqrt(dt);

% Integrazione numerica
% metodo: Heun stocastico
for n = 1:N

    % Costruzione del sistema al passo n-esimo
    M_n = M_S( S(:,n), gamma, alpha );  % M(S_n)
    drift  = M_n * H;                   
    diffus = M_n;                       

    % Generazione del rumore termico: processo di Wiener 3-d
    dW = sigma * sqrt_dt * randn(3,1);

    % Step predittore
    S_tilde = S(:,n) + drift * dt + diffus * dW;

    % Step correttore
    M_tilde = M_S( S_tilde, gamma, alpha ); % M(S_tilde)
    drift_tilde = M_tilde * H;
    diffus_tilde = M_tilde;

    % Step complessivo secondo Heun stocastico
    S(:,n+1) = S(:,n) + 0.5*(drift + drift_tilde)*dt + 0.5*(diffus + diffus_tilde)*dW;

    % Normalizzazione
    S(:,n+1) = S(:,n+1) / norm(S(:,n+1));

end

% Grafico componenti di spin
figure;
plot(t*1e12, [S(1,:); S(2,:); S(3,:)]','LineWidth', 2);
ax = gca; set(ax, 'FontSize', 14);
xlabel('Tempo [ps]', 'FontSize', 16);
ylabel('Componenti di spin', 'FontSize', 16);
legend('S_x','S_y','S_z', 'FontSize', 14, 'Location','best');
grid on; ax.GridAlpha = 0.3;

% Traiettoria 3D dello spin
figure;
plot3(S(1,:), S(2,:), S(3,:),'LineWidth', 2);
ax = gca; set(ax, 'FontSize', 14);
xlabel('S_x', 'FontSize', 16);
ylabel('S_y', 'FontSize', 16);
zlabel('S_z', 'FontSize', 16);
grid on; axis equal; view(135,30); % linea di vista dell'osservatore