clc
clear
close all

% Here we implement the pseudospectral method for the nonlinear
% Klein-Gordon Equation: u_tt -4u_xx + f(u) = 0, where f(u) = 4u -2u^3

% params:
omega = 1.98;
v= 0.1;
A = 2/(sqrt(3))*(sqrt(4-(omega)^2/(1-v^2)));
b = 2/sqrt(3)*(sqrt(1-v^2)/A);
k = (v*omega)/(1-v^2);

% now we set up our spatial grid: 
L = 100;
N = 512; % we recall that using a power of 2 for fft is most efficient. 
h = (2*L)/N; % since we are working on a spatial grid of [-L,L]...thus total length is 2L
x = (-N/2: N/2-1)*h;

% we find our Fourier wave numbers, and recall that we are working on an
% interval of [-L,L], thus we have 2pi.n/2L = (pi/L).n
q = (pi/L)*[0:N/2-1 -N/2:-1];

% we then setup our temporal grid: 
tau = 0.01;
T = 100;
Num_time_steps = round(T/tau);

% we introduce our initial conditions as follows: 
u_0 = A*cos(k*x).*sech(x/b);

u_t0 = A*(omega+(k*v))*sin(k*x).*sech(x/b)...
    + A*v/b*cos(k*x).*sech(x/b).*tanh(x/b);

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
U = [];
T_vals = [];

for n = 1:Num_time_steps
   
    u_hat = fft(u);
    u_xx = real(ifft(-(q.^2).*u_hat));
    
    % we use the centre difference to update in time: 
    u_new = 2*u - u_old + tau^2*(u_xx - 4*u + 2*u.^3);
    
    u_old = u;
    u = u_new;
    
    % we store the values within arrays in order to plot: 
    U = [U; u];
    T_vals = [T_vals; n*tau];
    
end
 
% finally we plot the solution: 
[X, TT] = meshgrid(x, T_vals);

figure;
surf(X, TT, U, 'EdgeColor', 'none',...
    'FaceColor', 'interp');
xlabel('x');
ylabel('t');
zlabel('u(x,t)');
title('Breather Solution Simulation using the Pseudospectral Method');
view(45,35);
colormap(turbo);
colorbar;
shading interp;
camlight;
lighting gouraud;