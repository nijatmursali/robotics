function q_dd = NS_acc(jac, sub, var3, var4, var5, var6, var7)
%02_KinematicRedundancy.pdf, page 58*
%Simply supply the arguments according to the formula, with composite
%terms decoupled first before going to the next term
%terms are added in the order 
%jac, sub, (r_dd, q_dot) | x_dd, q0_dd, K_D
%sub is an n * 2 matrix. First column holds the second argument of the
%subs matlab function, and the second column contains the third argument
%subs(jac, sub(:,1), sub(:,2))
%sub = [a  5]
%      [q1 p1/2]
%      [...     ]
%only the first four arguments are required, in which case, it defaults to a mode
%var6, if supplied, can only be a mode name; "mode1" or "mode2" 

%q_dd = NS_acc(jac, sub, x_dd, q0_dd)
%q_dd = NS_acc(jac, sub, r_dd, q_d, q0_dd)
%q_dd = NS_acc(jac, sub, x_dd, dqH, K_D, q_d, "mode1")
%q_dd = NS_acc(jac, sub, r_dd, q_d, dqH, K_D, "mode2")

% Robotics2Tools: https://github.com/MoSamArafat/Robotics2Tools
% Scripts based on Robotics 2 course at DIAG, Sapienza University of Rome taught by Prof. Alessandro De Luca.
% They are built to cover as many cases as found on previous exams.

% Most of the scripts have been tested on at least one example from previous exams.
% You are encouraged to test and provide feedback in case of any issues
 
% Contributors:
% Mohamed Samir Arafat
% Sayo Makinwa

J_d = pardiff(jac);

jac = subs(jac, sub(:,1), double(sub(:,2)));
jac = double(jac);

J_d = subs(J_d, sub(:,1), sub(:,2));
J_d = double(J_d);

[~, n] = size(jac);
J_inv = pinv(jac);

if (nargin == 4)
    x_dd = var3;
    q0_dd = var4;
    q_dd = (J_inv * x_dd) + ( (eye(n,n) - (J_inv * jac)) * q0_dd);
elseif (nargin == 5)
    r_dd = var3;
    q_d = var4;
    q0_dd = var5;
    x_dd = r_dd - (J_d * q_d);
    q_dd = (J_inv * x_dd) + ( (eye(n,n) - (J_inv * jac)) * q0_dd);
elseif (nargin == 7)
    if (var7 == "mode1")
        x_dd = var3;
        dqH = var4;
        K_D = var5;
        q_d = var6;
        q0_dd = dqH - (K_D * q_d);
        q_dd = (J_inv * x_dd) + ( (eye(n,n) - (J_inv * jac)) * q0_dd);
    elseif (var7 == "mode2")
        r_dd = var3;
        q_d = var4;
        dqH = var5;
        K_D = var6;
        x_dd = r_dd - (J_d * q_d);
        q0_dd = dqH - (K_D * q_d);
        q_dd = (J_inv * x_dd) + ( (eye(n,n) - (J_inv * jac)) * q0_dd);
    else
    end
end
end