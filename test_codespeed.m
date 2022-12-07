robot=Mccpvd1Dof();
model=robot.modelpara;
x=[0;0];
u=[0;0;0];
for i=1:100
    tic;
    %mccp1TauA(obj,xx,uu);
    robot.dynamics(x,u);
    Mtoc(i) = toc;
end

for j=1:100
    tic;
    mccp1Dynamics_mex(model,x,u);
     Mextoc(j) = toc;
end

t1 = mean(Mtoc)
t2 = mean(Mextoc)
t2/t1