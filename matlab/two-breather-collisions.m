clc
clear
close all

% params:
omega = 1.98;
v_1 = -0.1;   
v_2 =  0.1;   
A_1 = 2/sqrt(3)*sqrt(4-omega^2/(1-v_1^2));
b_1 = 2/sqrt(3)*sqrt(1-v_1^2)/A_1;
k_1 = v_1*omega/(1-v_1^2);

A_2 = 2/sqrt(3)*sqrt(4-omega^2/(1-v_2^2));
b_2 = 2/sqrt(3)*sqrt(1-v_2^2)/A_2;
k_2 = v_2*omega/(1-v_2^2);

d = 35;

% now we set up our spatial grid:
L = 100;
N = 512; % we recall that using a power of 2 for fft is most efficient.
h = 2*L/N; % since we are working on a spatial grid of [-L,L]...thus total length is 2L
x = (-N/2:N/2-1)*h;

% we find our Fourier wave numbers, and recall that we are working on an
% interval of [-L,L], thus we have 2pi.n/2L = (pi/L).n
q = (pi/L)*[0:N/2-1 -N/2:-1];

% we then setup our temporal grid: 
tau = 0.01;
T = 600;
Num_time_steps = round(T/tau);


% we introduce our initial conditions as follows: 
u_0 = A_1*cos(k_1*(x-d)).*sech((x-d)/b_1)...
   + A_2*cos(k_2*(x+d)).*sech((x+d)/b_2);

u_t0 = A_1*(omega+(k_1*v_1)).*sin(k_1*(x-d)).*sech((x-d)/b_1)...
    + A_1*(v_1/b_1).*cos(k_1*(x-d)).*sech((x-d)/b_1).*tanh((x-d)/b_1)...
    + A_2*(omega + k_2*v_2).*sin(k_2*(x+d)).*sech((x+d)/b_2)...
    + A_2*(v_2/b_2).*cos(k_2*(x+d)).*sech((x+d)/b_2).*tanh((x+d)/b_2);

% we then compute our spatial values in Fourier space, using u_0:
u_hat = fft(u_0);
u_xx0 = real(ifft(-(q.^2).*u_hat));

% next we compute u_tt(0) by subbing into our original equation: 
u_tt0 = u_xx0 - 4*u_0 + 2*u_0.^3;

% then we find the Taylor expansion of u(x,tau), which shall give us the
% next value of u in time, which we need for the central difference
% approximation, since the scheme is a 3 level time scheme: 
u_old = u_0;
u = u_0 + tau*u_t0 + (tau^2/2)*u_tt0;



% we instantiate a matrix U, which shall store the vectors of u for every
% point in time, and T_vals stores the corresponding time value for the
% time index of u: 
t_start_plot = 200;
U = [];
T_vals = [];

for n = 1:Num_time_steps

    u_hat = fft(u);
    u_xx = real(ifft(-(q.^2).*u_hat));

    % we use the centre difference to update in time: 
    u_new = 2*u - u_old + tau^2*(u_xx - 4*u + 2*u.^3);

    % we store the values within arrays in order to plot: 
    u_old = u;
    u = u_new;
    t=n*tau;

    % we truncate the starting point of some of our simulations for
    % computational efficiency: 
    if t >= t_start_plot
        U = [U; u];
        T_vals = [T_vals; t];
    end
end

% finally we plot the solution: 
[X, TT] = meshgrid(x, T_vals);

figure;
surf(X, TT, U, ...
    'EdgeColor', 'none', ...
    'FaceColor', 'interp');
xlabel('x');
ylabel('t');
zlabel('u(x,t)');
title('Collision of Two Breather Solutions with d = 50');
view(45,35);
colormap(turbo);
colorbar;
shading interp;
camlight;
lighting gouraud;

