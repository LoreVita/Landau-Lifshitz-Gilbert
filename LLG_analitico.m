%--------------------------------------------------------------------------
% Dinamica esatta dell'equazione di Landau-Lifshitz-Gilbert (LLG)
% Obiettivo: valutare la traiettoria di uno spin atomico classico 
%--------------------------------------------------------------------------

% Definizione della funzione sech
sech = @(x) 1 ./ cosh(x);

% Parametri
alpha = 0.1;     % smorzamento (adimensionale)
gamma = 1.76e11; % rapporto giromagnetico [rad/(s·T)]
H = 10;          % campo magnetico [T]

% Costanti
denom = 1 + alpha^2;
lambda = (alpha * gamma * H) / denom;
omega = (gamma * H) / denom;

% Intervallo temporale 
dt = 1e-15;           
t_max = 50e-12;      
t = 0:dt:t_max;     % vettore dei tempi

% Calcolo delle componenti di spin
Sx = sech(lambda * t) .* cos(omega * t);
Sy = sech(lambda * t) .* sin(omega * t);
Sz = tanh(lambda * t);

% Grafico componenti di spin
figure;
plot(t*1e12, [Sx; Sy; Sz], 'LineWidth', 2);
ax = gca; set(ax, 'FontSize', 14);
xlabel('Tempo [ps]', 'FontSize', 16);
ylabel('Componenti di spin', 'FontSize', 16);
legend('S_x(t)', 'S_y(t)', 'S_z(t)','FontSize', 14, 'Location','best');
grid on; ax.GridAlpha = 0.3;

% Traiettoria 3D dello spin
figure;
plot3(Sx, Sy, Sz, 'LineWidth', 2);
ax = gca; set(ax, 'FontSize', 14);
xlabel('S_x', 'FontSize', 16);
ylabel('S_y', 'FontSize', 16);
zlabel('S_z', 'FontSize', 16);
grid on; axis equal; view(135,30); % direzione di vista dell'osservatore