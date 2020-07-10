function q_dot = task_prio2 (jac, r_dot, map)
%02_KinematicRedundancy.pdf, page 51 - 56*; Robotics2_13.01.09.pdf, page 3

%jac is the jacobians for each task, stacked on top of one another,
%starting from highest proirity to the lowest
%r_dot follows intuitively from jac
%map is a [k * 2] matrix, where k is the number of tasks, kth row column 1
%contains the start positions in the stacked jac and r_dot, kth row column 2
%contains the end positions in the stacked jac and r_dot for the kth task

%we feel this function is more intuitively correct than the first one

% Robotics2Tools: https://github.com/MoSamArafat/Robotics2Tools
% Scripts based on Robotics 2 course at DIAG, Sapienza University of Rome taught by Prof. Alessandro De Luca.
% They are built to cover as many cases as found on previous exams.

% Most of the scripts have been tested on at least one example from previous exams.
% You are encouraged to test and provide feedback in case of any issues
 
% Contributors:
% Mohamed Samir Arafat
% Sayo Makinwa

[~, n] = size(jac);

J1_pinv = pinv(jac(map(1,1):map(1,2), :));
qk_dot_prev = J1_pinv * r_dot(map(1,1):map(1,2), :);
p_k_prev = eye(n) - (J1_pinv * jac(map(1,1):map(1,2), :));

[m, ~] = size(map);
for k = 2 : m
    jac_k = jac(map(k,1):map(k,2), :);
    rk_dot = r_dot(map(k,1):map(k,2), :);
    [qk_dot_prev, p_k_prev] = kth_task_prio (jac_k, qk_dot_prev, p_k_prev, rk_dot);
end
q_dot = qk_dot_prev;
end