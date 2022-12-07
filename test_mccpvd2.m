robot = Mccpvd2Dof();

x0 = zeros(4,1);
%u0 = [pi/3;pi/6;0.5;pi/3;pi/6;0.5];
u0 = [pi/2;pi/2;0.5;pi/2;pi/2;0.5];
%robot.step(x0,ones(6,1),0.02)
u = repmat(u0,1,100);

xsim = robot.simulate_feedforward(x0,u,0.02);

%%
%global mccpvd2_fh mccpvd2_ph1 mccpvd2_ph2
robot.plot(x0);

%%
robot.update(xsim(:,end));
%%

robot.plot([-pi/3; 2*pi/3])

%%
frames = robot.animate(xsim,0.02);
robot.savegif(frames, 'mccpvd2animation')
%%
movie(frames,1,50)

%%
gifname = 'mccpvd2_animation.gif';
for idx = 1:length(frames)
    im = frame2im(frames(idx));
    [A,map] = rgb2ind(im,256);
    if idx == 1
        imwrite(A,map,gifname,'gif','LoopCount',Inf,'DelayTime',0.02);
    else
        imwrite(A,map,gifname,'gif','WriteMode','append','DelayTime',0.02);
    end
end

