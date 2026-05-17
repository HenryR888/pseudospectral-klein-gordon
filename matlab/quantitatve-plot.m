clc
clear
close all

% params:
omega = 1.98;
v_1 = -0.1;   
v_2 =  0.1;   

% spatial grid:
L = 100;
N = 512;
h = 2*L/N;
x = (-N/2:N/2-1)*h;
q = (pi/L)*[0:N/2-1 -N/2:-1];
tau = 0.01;

% we test the following values of d: 
d_values = 0:5:50;

% here we create the arrays to store our amplitudes for collision and
% times:
collision_amplitudes = zeros(size(d_values));
collision_times = zeros(size(d_values));

for idx = 1:length(d_values)

    d = d_values(idx);

    % params for first breather: 
    A_1 = 2/sqrt(3)*sqrt(4-omega^2/(1-v_1^2));
    b_1 = 2/sqrt(3)*sqrt(1-v_1^2)/A_1;
    k_1 = v_1*omega/(1-v_1^2);

    % params for the second breather: 
    A_2 = 2/sqrt(3)*sqrt(4 - omega^2/(1-v_2^2));
    b_2 = 2/sqrt(3)*sqrt(1-v_2^2)/A_2;
    k_2 = v_2*omega/(1-v_2^2);

    T_collision_estimate = d/abs(v_2);

    if d == 0
        T = 100;
        T_collision_estimate = 0;
    else
        T = T_collision_estimate + 100;
    end

    Num_time_steps = round(T/tau);

    % ICs: 
    u_0 = A_1*cos(k_1*(x-d)).*sech((x-d)/b_1) ...
        + A_2*cos(k_2*(x+d)).*sech((x+d)/b_2);

    u_t0 = A_1*(omega + k_1*v_1).*sin(k_1*(x-d)).*sech((x-d)/b_1) ...
        + A_1*(v_1/b_1).*cos(k_1*(x-d)).*sech((x-d)/b_1).*tanh((x-d)/b_1) ...
        + A_2*(omega + k_2*v_2).*sin(k_2*(x+d)).*sech((x+d)/b_2) ...
        + A_2*(v_2/b_2).*cos(k_2*(x+d)).*sech((x+d)/b_2).*tanh((x+d)/b_2);

    u_hat = fft(u_0);
    u_xx0 = real(ifft(-(q.^2).*u_hat));

    u_tt0 = u_xx0 - 4*u_0 + 2*u_0.^3;

    u_old = u_0;
    u = u_0 + tau*u_t0 + (tau^2/2)*u_tt0;

    max_collision_amp = 0;
    time_of_max_amp = 0;

    collision_window = 30;

    for n = 1:Num_time_steps

        t = n*tau;

        u_hat = fft(u);
        u_xx = real(ifft(-(q.^2).*u_hat));

        u_new = 2*u - u_old + tau^2*(u_xx - 4*u + 2*u.^3);

        u_old = u;
        u = u_new;

        if abs(t - T_collision_estimate) <= collision_window

            current_amp = max(abs(u));

            if current_amp > max_collision_amp
                max_collision_amp = current_amp;
                time_of_max_amp = t;
            end

        end

    end

    collision_amplitudes(idx) = max_collision_amp;
    collision_times(idx) = time_of_max_amp;

    fprintf('d = %.1f, collision time approx = %.2f, amplitude = %.5f\n', ...
        d, time_of_max_amp, max_collision_amp);

end

figure;
plot(d_values, collision_amplitudes, 'o-', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('Half-separation distance d');
ylabel('Amplitude at collision');
title('Collision Amplitude Against Half-Separation Distance d');
grid on;