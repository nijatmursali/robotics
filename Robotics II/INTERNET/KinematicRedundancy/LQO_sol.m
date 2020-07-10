function [q_dot, pr_mx] = LQO_sol(jac, q0_dot, W, r_dot)
%02_KinematicRedundancy.pdf, page 29 - 31 - 32*
%conversion to rad from deg: deg*pi/180
% W must be a square matrix [nJ x nJ]

% Robotics2Tools: https://github.com/MoSamArafat/Robotics2Tools
% Scripts based on Robotics 2 course at DIAG, Sapienza University of Rome taught by Prof. Alessandro De Luca.
% They are built to cover as many cases as found on previous exams.

% Most of the scripts have been tested on at least one example from previous exams.
% You are encouraged to test and provide feedback in case of any issues
 
% Contributors:
% Mohamed Samir Arafat
% Sayo Makinwa

jac = double(jac);

W_inv = inv(W);
Jw_a = W_inv * jac';
Jw_b = inv(jac * W_inv * jac');
Jw = Jw_a * Jw_b;
%q_dot = q_dot0 - (Jw * jac*q_dot0 - r_dot );
[m, ~] = size(Jw);
pr_mx = eye(m) - (Jw * jac);
q_dot = Jw * r_dot + (pr_mx*q0_dot);

end