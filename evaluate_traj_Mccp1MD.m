function [ traj ] = evaluate_traj_Mccp1MD( robot, traj, varargin )
% evaluate energy and performance of

parser = inputParser();
addOptional(parser, 'qf', NaN);
parse(parser, varargin{:});


traj.stiffness = robot.stiffness( traj.x );
traj.damping = robot.damping(traj.u);
sdt = 0.001;
tsim = 0:sdt:traj.t(end); % simulation timestamp
usim = scale_controlSeq(traj.u, traj.t, tsim);
psim.solver = 'rk4';
psim.dt = sdt;
f = @(x,u) robot.dynamics(x,u);
xsim = simulate_feedforward(f, traj.x(:,1), usim, psim);
traj.tsim = tsim;
traj.usim = usim;
traj.xsim = xsim;
traj.qf = parser.Results.qf;
if isnan(traj.qf)
    sinfo = stepinfo(traj.xsim(1,:), traj.tsim, 'SettlingTimeShreshold', 0.02);
else
    sinfo = stepinfo(traj.xsim(1,:), traj.tsim, traj.qf, 'SettlingTimeShreshold', 0.02);
end
traj.SettlingTime = sinfo.SettlingTime;
traj.Pout = robot.power_out(xsim(:,1:end-1),usim);
traj.Pinp_ms = robot.power_in(xsim(:,1:end-1), usim);
traj.Pin1 = robot.power_in1(xsim(:,1:end-1), usim);
traj.Pin1p = max(traj.Pin1,0);
traj.Pin2 = robot.power_in2(xsim(:,1:end-1), usim);
traj.Pin2p = max(traj.Pin2,0);
traj.power_rege = robot.power_rege(xsim(:,1:end-1),usim);
traj.power_damp = robot.power_damp(xsim(:,1:end-1),usim);
[traj.Pe, traj.Pe1, traj.Pe2] = robot.power_elec(xsim(:,1:end-1),usim, 'MotoringOnly', 1);
%traj.sqr_torque = robot.sqrt_torques(xsim(:,1:end-1),usim);
traj.Es = robot.energy_spring(xsim,usim);
traj.Esd = gradient(traj.Es,sdt);
traj.Pinp = max(traj.Esd(1:end-1) + traj.Pout,0);
traj.Pin = traj.Esd(1:end-1) + traj.Pout;

traj.Eout = sum(traj.Pout, 2)*sdt;
traj.Ein = sum(traj.Pin, 2)*sdt;
traj.Einp = sum(traj.Pinp, 2)*sdt;
traj.Einp_ms = sum(traj.Pinp_ms, 2)*sdt;
%traj.Ein_net = sum( traj.Esd(1:end-1) + traj.Pout, 2)*sdt;
%traj.Ein1 = sum(power_in1p)*sdt;
%traj.Ein2 = sum(power_in2p)*sdt;
traj.Edamp = sum(traj.power_damp,2)*sdt;
traj.Erege = sum(traj.power_rege,2)*sdt;
traj.Ee = sum(traj.Pe)*sdt;
traj.Ee1 = sum(traj.Pe1)*sdt;
traj.Ee2 = sum(traj.Pe2)*sdt;
%traj.Esqr = sum(sqr_torque)*sdt;

end

