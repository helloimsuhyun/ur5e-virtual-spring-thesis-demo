function tau = ur5e__ctrl(u)
% u = [ q qd x_ref ] 6 6 3  → 15

persistent robot eeName R_ref_fixed Kp Kd K_R K_w p_home

if isempty(robot)
    robot  = loadrobot('universalUR5e','DataFormat','row','Gravity',[0 0 -9.81]);
    eeName = 'tool0';

    q_home = [0 -pi/2 pi/2 0 pi/2 0];
    T0 = getTransform(robot, q_home, eeName);
    p_home = T0(1:3,4);
    R_ref_fixed = [ 0 0 1;
                    0 1 0;
                   -1 0 0 ];
    % gain
    % 위치 제어 (position)
    Kp = diag([360 360 360]);    
    Kd = diag([120 120 120]);      

    % 자세 제어 (orientation)
    K_R = diag([80 80 80]);     
    K_w = diag([20 20 20]);   


end

u = u(:).';        % 1x15
q   = u(1:6);
qd  = u(7:12);
x_ref = u(13:15).';    % 3x1

% FK
T = getTransform(robot, q, eeName);
p = T(1:3,4);          % 3x1
R = T(1:3,1:3);        % 3x3

% Jacobian
J  = geometricJacobian(robot, q, eeName);
Jw = J(1:3,:);         % 3x6
Jv = J(4:6,:);         % 3x6 geometricJacobian은 위에 3개가 Jw임

%% task-space- virtual spring-damper
x_dot = Jv * qd(:);               % 3x1
e_x   = x_ref - p;                % 3x1
F_pos = Kp * e_x - Kd * x_dot;    % 3x1

%% --- orientation error R  ---
xd = R_ref_fixed(:,1); yd = R_ref_fixed(:,2); zd = R_ref_fixed(:,3);
x  = R(:,1);           y  = R(:,2);           z  = R(:,3);

e_R = 0.5 * ( cross(x,xd) + cross(y,yd) + cross(z,zd) );  % orientation error 3x1

omega = Jw * qd(:);              % 현재 각속도
omega_ref = zeros(3,1);          % 목표 각속도 0 → 자세 고정
e_w = omega_ref - omega;

M_rot = K_R * e_R + K_w * e_w;   % 3x1
%M_rot = zeros(3,1);   % 자세 제어 비활성화


F_task = [M_rot ; F_pos];         % 6x1
tau_t = (J.' * F_task).';
tau_joint_damp = -1.0 * qd;  % joint space damping term

% gravity compensation term
tau_g = gravityTorque(robot, q);  % 1x6

%% --- control input u
tau = tau_t + tau_g + tau_joint_damp;  % 1x6

end
