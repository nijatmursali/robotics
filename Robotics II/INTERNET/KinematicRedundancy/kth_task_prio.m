function [qk_dot, p_k] = kth_task_prio (jac_k, qk_dot_prev, p_k_prev, rk_dot)
%02_KinematicRedundancy.pdf, page 56*

% Robotics2Tools: https://github.com/MoSamArafat/Robotics2Tools
% Scripts based on Robotics 2 course at DIAG, Sapienza University of Rome taught by Prof. Alessandro De Luca.
% They are built to cover as many cases as found on previous exams.

% Most of the scripts have been tested on at least one example from previous exams.
% You are encouraged to test and provide feedback in case of any issues
 
% Contributors:
% Mohamed Samir Arafat
% Sayo Makinwa

Jk_pinv = pinv(jac_k);
[~, n] = size(jac_k);

p_k = eye(n) - (Jk_pinv * jac_k);

qk_dot = qk_dot_prev + ( pinv(jac_k * p_k_prev) * (rk_dot - (jac_k * qk_dot_prev) ));
disp('i been called');
end