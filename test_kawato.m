%clear all

% set up simulation
dt = .02;              % control frequency
N  = 100;              % no. samples
t  = (0:N-1)*dt;       % sample times
x0 = [pi/3;pi/3;0;0];  % start state

% set up dynamics function
%m  = model_kawato_arm; % model parameters
arm = Kawato2Dof();
%f = @(x, u) arm.dynamics  ( x, u); % dynamics function

% simulation parameter struct
%p = [];
dt = 0.02;
%p.solver = 'rk4';

% define feed forward commands
U = sin(5*rand(6,1)*t);

% simulate
tic
X = arm.simulate_feedforward ( x0, U, dt );
%运动学部分
%X(1) = theta1  X(2) = theta2
%求v，h 到x


toc
%%
% plot joint angles, velocities, commands
figure
subplot(3,1,1),plot(X(1:2,:)'),ylabel('q (rad)')
ylim([-pi/2 pi/2])
subplot(3,1,2),plot(X(3:4,:)'),ylabel('dq/dt (rad/s)')
subplot(3,1,3),plot(U'),ylabel('u'),xlabel('t (s)')

arm.animate(X,0.02);