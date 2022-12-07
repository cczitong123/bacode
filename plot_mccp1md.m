function [ figh ] = plot_mccp1md( traj, varargin )
%PLOT_MCCP1 Summary of this function goes here
%   Detailed explanation goes here

parser = inputParser();
%addRequired(parser, 'data');
addOptional(parser, 'PlotTarget', 0);
addOptional(parser, 'PlotSettlingTime', 0);
addOptional(parser, 'PlotEnergy', 0);
addOptional(parser, 'PlotImpedance', 0);
addOptional(parser, 'FigHandle', []);
parse(parser, varargin{:});
opts = parser.Results;

t = traj.t;
if isempty(opts.FigHandle)
figh = figure('Position', [100 100 550 950]);
else
    figh = figure(opts.FigHandle);
end
subplot(311)
hold on
plot(t, traj.x(1,:))
plot(t, traj.x(3,:),'--')
title('Joint & EP trajectories')
legend('Joint','EP')
hold off
subplot(312)
hold on
plot(t, traj.x(4,:))
ylim([0 pi/2])
if isfield(traj, 'stiffness')
    yyaxis right
    plot(t(1:length(traj.stiffness)), traj.stiffness,'--')
    title('Stiffness Motor trajectory & Stiffness profile')
else
    title('Stiffness Motor trajectory')
end

hold off
subplot(313)
hold on
plot(t(1:end-1), traj.u(3,:))
title('Damping')
hold off

if opts.PlotEnergy == 1
figure
subplot(211)
hold on
plot( traj.t(1:length(traj.Pin)), traj.Pin )
plot( traj.t(1:length(traj.Pout)), traj.Pout )
plot( traj.t(1:length(traj.Esd)), traj.Esd )
plot( traj.t(1:length(traj.Pin_net)), traj.Pinp )
hold off
legend('Pin','Pout','Esd','[Pin]+')

subplot(212)
plot( traj.t(1:length(traj.Pin1)),traj.Pin1)
plot( traj.t(1:length(traj.Pin2)),traj.Pin2)
plot( traj.t(1:length(traj.Pinp_ms)),traj.Pinp_ms)
end

if opts.PlotImpedance == 1
    figure
    subplot(211)
    plot(traj.t(1:length(traj.stiffness)), traj.stiffness)
    title('Stiffness')
    ylabel('K (Nm/rad)')
    subplot(212)
    plot(traj.t(1:length(traj.damping)), traj.damping)
    title('Damping')
    ylabel('D (Nms/rad)')
    xlabel('Time (s)')
end

end