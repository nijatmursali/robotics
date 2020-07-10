function q_dSNS = SNS_sol(jac, x_dot, q_lim)
%function a = SNS_sol(jac, x_dot, q_lim)
%     clc
%02_KinematicRedundancy.pdf, page 58*

%q_dd = NS_acc(jac, sub, x_dd, q0_dd)

% Robotics2Tools: https://github.com/MoSamArafat/Robotics2Tools
% Scripts based on Robotics 2 course at DIAG, Sapienza University of Rome taught by Prof. Alessandro De Luca.
% They are built to cover as many cases as found on previous exams.

% Most of the scripts have been tested on at least one example from previous exams.
% You are encouraged to test and provide feedback in case of any issues
 
% Contributors:
% Mohamed Samir Arafat
% Sayo Makinwa

jac = double(jac);
J_inv = pinv(jac);
q_dot = J_inv * x_dot;
q_len = length(q_lim);

W = eye(q_len);
q_dN = zeros(q_len,1);
s = 0;

[mj, nj] = size(jac);
for (k = 1: (nj-mj))

    q_dot = ((pinv(jac*W)) * (1 * x_dot)) + ((eye(q_len) - (pinv(jac*W)*jac)) * q_dN)
    q_temp = q_lim - abs(q_dot);
    temp_index = find(q_temp == unique(min(q_temp)));
    s = max((q_lim(temp_index(1))/q_dot(temp_index(1))), s)

    saturated = [];
    for i = 1:q_len
        if (q_lim(i) <= abs(q_dot(i)))
            saturated = [saturated; [i, q_dot(i)/abs(q_dot(i))]]
        end
    end

    if isempty(saturated)
        break;
    end

    [m,n] = size(saturated);
    for i = 1:m
        W(saturated(i,1), saturated(i,1)) = 0;
        q_dN(saturated(i,1)) = q_lim(saturated(i,1)) * saturated(i,2);
    end
end

pre_q_dSNS = ((pinv(jac*W)) * (1 * x_dot)) + ((eye(q_len) - (pinv(jac*W)*jac)) * q_dN)

if (min(q_temp) >= 0)
    q_dSNS = pre_q_dSNS;
else
    q_temp = q_lim - abs(pre_q_dSNS);
    temp_index = find(q_temp == min(q_temp));
    s = max((q_lim(temp_index(1))/q_dot(temp_index(1))), s)
    q_dSNS = ((pinv(jac*W)) * (s * x_dot)) + ((eye(q_len) - (pinv(jac*W)*jac)) * q_dN);
end
end