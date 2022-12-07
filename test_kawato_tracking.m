
x0=[0;0];
%y0=0;
xf = [pi/3;pi/3];
%yf = pi/3;

[x, t] = generate_trajectory_jerk(x0,xf,0.7,0.02);
%y = generate_trajectory_jerk(y0,yf,0.7,0.02);

