function gingerWaistIKConstrained()
close all;
% clear; clc;
rosshutdown;

targetPos = [0.2; -0.2; 0.25]; % [1:3] Postion, [4:6] AngleAxis
targetAA = [0, 0, 1, 30/180*pi];

rosinit;
JSpub = rospublisher("/joint_states","sensor_msgs/JointState"); 
JSmsg = rosmessage(JSpub);
JSmsg.Name={'Back_Z', 'Back_X', 'Back_Y', 'Left_Shoulder_X', 'Left_Shoulder_Y', 'Left_Elbow_Z', 'Left_Elbow_X', 'Left_Wrist_Z', 'Left_Wrist_X', 'Left_Wrist_Y'};

signM = ones(1,10)*(-1);
signM(1) = 1; signM(3) = 1; signM(5) = 1; signM(10) = 1;
%%
th0 = 20.0/180.0*pi;
c =[
    % waist
    -1.54, 1.54;
    -0.5, 0.2;
    -0.38, 0.38;
    % shoudlder
    -0.6, 3.14;
    -0.38, 1.54;
    -1.54, 1.54;
    % elbow
    -0.2, 2.1;
    % wrist
    -1.54, 1.54;
    -0.38, 0.38;
    -0.6, 0.38;
]; %constraint



%%
syms uth1 uth2 uth3 uth4 uth5 uth6 uth7 uth8 uth9 uth10; %under theta
utheta = [uth1; uth2; uth3; uth4; uth5; uth6; uth7; uth8; uth9; uth10];

th1_ =  (c(1,2)-c(1,1))/pi*atan(uth1)+ (c(1,1)+c(1,2))/2;
th2_ =  (c(2,2)-c(2,1))/pi*atan(uth2)+ (c(2,1)+c(2,2))/2;
th3_ =  (c(3,2)-c(3,1))/pi*atan(uth3)+ (c(3,1)+c(3,2))/2;
th4_ =  (c(4,2)-c(4,1))/pi*atan(uth4)+ (c(4,1)+c(4,2))/2;
th5_ =  (c(5,2)-c(5,1))/pi*atan(uth5)+ (c(5,1)+c(5,2))/2;
th6_ =  (c(6,2)-c(6,1))/pi*atan(uth6)+ (c(6,1)+c(6,2))/2;
th7_ =  (c(7,2)-c(7,1))/pi*atan(uth7)+ (c(7,1)+c(7,2))/2;
th8_ =  (c(8,2)-c(8,1))/pi*atan(uth8)+ (c(8,1)+c(8,2))/2;
th9_ =  (c(9,2)-c(9,1))/pi*atan(uth9)+ (c(9,1)+c(9,2))/2;
th10_ = (c(10,2)-c(10,1))/pi*atan(uth10)+ (c(10,1)+c(10,2))/2;
theta = [th1_; th2_; th3_; th4_; th5_; th6_; th7_; th8_; th9_; th10_];
J_ =jacobian(theta, utheta);

syms th1 th2 th3 th4 th5 th6 th7 th8 th9 th10;
theta = [th1; th2; th3; th4; th5; th6; th7; th8; th9; th10];

% waist
M01 = [1, 0, 0, 0;
    0, cos(th1), -sin(th1), 0;
    0, sin(th1), cos(th1), 0;
    0, 0, 0, 1;];
M12 = [cos(th2), 0, sin(th2), 0;
    0, 1, 0, 0;
    -sin(th2), 0, cos(th2), 0;
    0, 0, 0, 1;];
M23 = [cos(th3), -sin(th3), 0, 0.105;
    sin(th3), cos(th3), 0, 0;
    0, 0, 1, 0;
    0, 0, 0, 1; ] * ...
    [cos(th0), -sin(th0), 0, 0.25743421;
    sin(th0), cos(th0), 0, -0.167643;
    0, 0, 1, -0.01532;
    0, 0, 0, 1 ];

% shoudler
M34 =  [cos(th4), 0,  sin(th4), 0;
    0, 1, 0, 0;
    -sin(th4), 0, cos(th4), 0;
    0, 0, 0, 1 ];
M45 = [cos(th5), -sin(th5), 0, 0;
    sin(th5), cos(th5), 0, 0;
    0, 0, 1, 0;
    0, 0, 0, 1 ] * ...
    [cos(th0), -sin(-th0), 0, -0.22791;
    sin(-th0), cos(th0), 0, 0.03067462;
    0, 0, 1, 0;
    0, 0, 0, 1 ] ;
M56 = [1, 0, 0, 0;
    0, cos(th6), -sin(th6), 0;
    0, sin(th6), cos(th6), 0;
    0, 0, 0, 1];

% elbow
M67 = [cos(th7), 0, sin(th7),  0;
    0, 1, 0, 0;
    -sin(th7), 0, cos(th7),  0;
    0, 0, 0, 1 ] * ...
    [1, 0, 0, -0.209557780;
    0, 1, 0, 0;
    0, 0, 1, 0;
    0, 0, 0, 1];

% wrist
M78 = [1, 0, 0, 0;
    0, cos(th8), -sin(th8), 0;
    0, sin(th8), cos(th8), 0;
    0, 0, 0, 1 ];
M89 = [cos(th9), 0,  -sin(th9), 0;
    0, 1, 0, 0;
    sin(th9), 0, cos(th9), 0;
    0, 0, 0, 1 ];
M910 = [cos(th10), -sin(th10), 0, 0;
    sin(th10), cos(th10), 0, 0;
    0, 0, 1, 0;
    0, 0, 0, 1 ];

M010 = M01*M12*M23*M34*M45*M56*M67*M78*M89*M910;

M10_9 = inv(M910);
M10_8 = M10_9 / M89;
M10_7 = M10_8 / M78;
M10_6 = M10_7 / M67;
M10_5 = M10_6 / M56;
M10_4 = M10_5 / M45;
M10_3 = M10_4 / M34;
M10_2 = M10_3 / M23;
M10_1 = M10_2 / M12;
M10_0 = M10_1 / M01;

% pos = M010(1:3, 4);
% aa = rotm2axang(M010(1:3,1:3));
% Jpos = jacobian(pos, theta);
% Jpos = vpa(Jpos, 5);


%xyz yzx y xyz
Jfull_10 =  [
    [cross(M10_0(1:3, 4), M10_0(1:3, 1)); M10_0(1:3, 1)], ...
    [cross(M10_1(1:3, 4), M10_1(1:3, 2)); M10_1(1:3, 2)], ...
    [cross(M10_2(1:3, 4), M10_2(1:3, 3)); M10_2(1:3, 3)], ...
    [cross(M10_3(1:3, 4), M10_3(1:3, 2)); M10_3(1:3, 2)], ...
    [cross(M10_4(1:3, 4), M10_4(1:3, 3)); M10_4(1:3, 3)], ...
    [cross(M10_5(1:3, 4), M10_5(1:3, 1)); M10_5(1:3, 1)], ...
    [cross(M10_6(1:3, 4), M10_6(1:3, 2)); M10_6(1:3, 2)], ...
    [cross(M10_7(1:3, 4), M10_7(1:3, 1)); M10_7(1:3, 1)], ...
    [cross(M10_8(1:3, 4), M10_8(1:3, 2)); M10_8(1:3, 2)], ...
    [cross(M10_9(1:3, 4), M10_9(1:3, 3)); M10_9(1:3, 3)]  ...
];

Jfull_0 = [M010(1:3, 1:3), zeros(3); zeros(3),M010(1:3, 1:3)] * Jfull_10;

% vpa(simplify(Jfull_0), 5)
% return 

%% 
% curUTheta = [
%     tan( ((c(1,1)+c(1,2))*rand()-(c(1,1)+c(1,2))/2 - (c(1,1)+c(1,2))/2) * pi/(c(1,2)-c(1,1)));
%     tan( ((c(2,1)+c(2,2))*rand()-(c(2,1)+c(2,2))/2 - (c(2,1)+c(2,2))/2) * pi/(c(2,2)-c(2,1)));
%     tan( ((c(3,1)+c(3,2))*rand()-(c(3,1)+c(3,2))/2 - (c(3,1)+c(3,2))/2) * pi/(c(3,2)-c(3,1)));
%     tan( ((c(4,1)+c(4,2))*rand()-(c(4,1)+c(4,2))/2 - (c(4,1)+c(4,2))/2) * pi/(c(4,2)-c(4,1)));
% ];

% curUTheta = [
%     tan( (0.1 - (c(1,1)+c(1,2))/2) * pi/(c(1,2)-c(1,1)));
%     tan( (0.1 - (c(2,1)+c(2,2))/2) * pi/(c(2,2)-c(2,1)));
%     tan( (0.1 - (c(3,1)+c(3,2))/2) * pi/(c(3,2)-c(3,1)));
%     tan( (0.1 - (c(4,1)+c(4,2))/2) * pi/(c(4,2)-c(4,1)));
% ];

%%
curTheta = ones(10, 1)*0.1;
% curTheta = zeros(10, 1);
% curTheta(1:3) = [0;0;0];

M010_value = double(subs(M010, theta, curTheta));
curPos = M010_value(1:3,4);
v = targetPos-curPos;
v = [0;0;0];
curM = M010_value(1:3,1:3);
deltaR = curM' * axang2rotm(targetAA) ;
deltaAA = rotm2axang(deltaR);
w = deltaAA(1:3)'*deltaAA(4);
V = [v;  w];

dt = 0.5; 
idx = 0;

% while norm(target-curPose) > 0.01 

while norm(w) > 0.01
    idx = idx+1;
%     X = sprintf('%d: norm: %f  %f.', idx, norm(v), norm(w));
%     disp(X);
    
    myJfull_0 = double(subs(Jfull_0, theta, curTheta));
    myJfull_10 = double(subs(Jfull_10, theta, curTheta));
    myInvJfull_0 = pinv(myJfull_0)
    
    deltaTheta = myInvJfull_0 * V * dt;
    curTheta = curTheta + deltaTheta;

    M010_value = double(subs(M010, theta, curTheta));
    curPos = M010_value(1:3,4);
    v = targetPos-curPos;
    v = [0;0;0];
    curM = M010_value(1:3,1:3);
    deltaR = curM' * axang2rotm(targetAA);
    deltaAA = rotm2axang(deltaR);
    w = deltaAA(1:3)'*deltaAA(4)/10;
    V = [v; w];
    
    angs = double(vpa(curTheta, 5))' .* signM;
    JSmsg.Position = angs;
    send(JSpub,JSmsg);
end

curM;
TargetM = axang2rotm(targetAA);

end