function [ xdot ] = mccp1Dynamics( model, x, u )
%MCCP1DYNAMICS Summary of this function goes here
%   Detailed explanation goes here
qdd = accel(model, x, u);
xdot = [x(2,:); qdd];

end

function qdd = accel(model, x, u)
tau = tauA(model, x, u);
qdd = (tau - model.fv*x(2,:) + gravity(x(1,:), model.mass, model.com, model.gc))/model.M;
end

function tau = tauA(model, x, u)
tau = model.Ks*model.B*model.C*sin(u(1,:)-x(1,:)).*(1+ (model.r*u(2,:) ...
                -model.A0)./sqrt(model.B^2+model.C^2-2*model.B*model.C*cos(u(1,:)-x(1,:))) )...
                - model.Dm*u(3,:).*x(2,:);
end

function f = gravity(q, mass, com, gc)
f = -mass*com*gc*sin(q);
end