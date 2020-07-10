function [q_dot, pr_mx] = PG_sol(jac, q0_dot, r_dot)
%02_KinematicRedundancy.pdf, page 33*
%conversion to rad from deg: deg*pi/180

% Robotics2Tools: https://github.com/MoSamArafat/Robotics2Tools
% Scripts based on Robotics 2 course at DIAG, Sapienza University of Rome taught by Prof. Alessandro De Luca.
% They are built to cover as many cases as found on previous exams.

% Most of the scripts have been tested on at least one example from previous exams.
% You are encouraged to test and provide feedback in case of any issues
 
% Contributors:
% Mohamed Samir Arafat
% Sayo Makinwa

jac = double(jac);

J_pinv = pinv(jac);
q_dot = q0_dot + (J_pinv * (r_dot - (jac * q0_dot)) );

[~, n] = size(jac);
pr_mx = eye(n) - (J_pinv * jac);

end