
f0 = 1000; %Nyquist en Hz 
alphas = [0, 0.25, 0.75, 1]; %roll-off indicado
t = linspace(0, 0.02, 1000); % Tiempo >= 0 (impulso)
f = linspace(-2*f0*(1+1), 2*f0*(1+1), 2000); % Frecuencia de -2B a 2B

% Figura respuestas al impulso
figure;
for i = 1:length(alphas)
    alpha = alphas(i);
    f_delta = alpha * f0;
    B = f0 + f_delta;
    % Calcular respuesta al impulso h_e(t)
    sinc_term = sin(2*pi*f0*t) ./ (2*pi*f0*t);
    sinc_term(t==0) = 1; % corrige división por 0
    cos_term = cos(2*pi*f_delta*t);
    denom = 1 - (4*f_delta*t).^2;
    impulse_response = 2*f0 * sinc_term .* (cos_term ./ denom);%funcion final
    subplot(length(alphas), 1, i)
    plot(t, impulse_response, 'LineWidth', 1.5);
    title(['Respuesta al impulso alfa = ' num2str(alpha)]);
    xlabel('Tiempo (s)');
    ylabel('h_e(t)');
    grid on;
end

% Figura respuestas en frecuencia
figure;
for i = 1:length(alphas)
    alpha = alphas(i);
    f_delta = alpha * f0;
    B = f0 + f_delta;
    f1 = f0 - f_delta;
    He_f = zeros(size(f));
    abs_f = abs(f);
    for j = 1:length(f)%funcion por tramos / funcion 10
        if abs_f(j) < f1
            He_f(j) = 1;
        elseif abs_f(j) < B
            He_f(j) = 0.5 * (1 + cos(pi*(abs_f(j) - f1)/(2*f_delta)));
        else
            He_f(j) = 0;
        end
    end
    subplot(length(alphas), 2, i)
    plot(f, He_f, 'LineWidth', 1.5);
    title(['Respuesta en frecuencia alfa = ' num2str(alpha)]);
    xlabel('Frecuencia (Hz)');
    ylabel('H_e(f)');
    grid on;
end

%espacio

% 1. Respuestas al impulso (coseno alzado)
f0 = 9000; % Frecuencia de Nyquist en Hz
alphas = [0, 0.25, 0.75, 1]; % Roll-off
t = linspace(-0.005, 0.005, 1000); % Tiempo centrado
impulse_responses = cell(length(alphas), 1);

for i = 1:length(alphas)
    alpha = alphas(i);
    f_delta = alpha * f0;
    sinc_term = sin(2*pi*f0*t) ./ (2*pi*f0*t); %formula 14
    sinc_term(t == 0) = 1; % evitar NaN
    cos_term = cos(2*pi*f_delta*t);
    denom = 1 - (4*f_delta*t).^2;
    impulse_response = 2*f0 * sinc_term .* (cos_term ./ denom);
    impulse_responses{i} = impulse_response;
end


% 2. Señal de entrada 

symbols = randi([0 1], 100, 1)*2 - 1; 
sps = 10; % muestras
upsampled = zeros(1, length(symbols) * sps);
upsampled(1:sps:end) = symbols;


% 3. Convolución con el filtro (alpha = 0.25)

selected_alpha_idx = 4 ; % Corresponde a alpha = 0.25
% Corresponde a 1 = 0, 2 = 0.25, 3 = 0.75, 4 = 1, se debe escoger que valor
% de a usar entre los 4
h = impulse_responses{selected_alpha_idx};
tx_signal = conv(upsampled, h, 'same');


% 4.ruido blanco AWGN

SNR = 18; 
tx_signal_noisy = awgn(tx_signal, SNR, 'measured');  


% 5. Diagrama de ojo

figure;
eyediagram(tx_signal_noisy, 2*sps); % 2 símbolos por trazo
ylim([-0,5 0,5]); % ajustar según amplitud
title('Diagrama de Ojo para \alpha = 0.25 con Ruido');