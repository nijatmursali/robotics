function q_dot = task_aug_sol (jac, jac_y, r_dot, y_dot, W_or_q0_dot)
% [gj, hj] = size (jac);
% [gy, hy] = size (jac_y);
% [gr ,hr] = size (r_dot);
% [gy, hy] = size (y_dot);
% [gwq, hwq] = size (W_or_q0_dot);

% Robotics2Tools: https://github.com/MoSamArafat/Robotics2Tools
% Scripts based on Robotics 2 course at DIAG, Sapienza University of Rome taught by Prof. Alessandro De Luca.
% They are built to cover as many cases as found on previous exams.

% Most of the scripts have been tested on at least one example from previous exams.
% You are encouraged to test and provide feedback in case of any issues
 
% Contributors:
% Mohamed Samir Arafat
% Sayo Makinwa

Ja = [jac; jac_y];
[~, n] = size(Ja);

% if (isequal([gr ,hr],[gj, 1])) %r_dot dimension check
%     ;
% elseif (isequal([gr ,hr], [1, gj]))
%     r_dot = r_dot';
% else
%     disp('Incorrect dimensions of either r_dot or jac')
% end
% 
% if (isequal([gwq ,hwq],[n, 1])) %q0_dot dimension check
%     q0_dot = W_or_q0_dot;
%     W_or_q0_dot = NaN;
% elseif (isequal([gwq ,hwq],[hj, 1]))
%     q0_dot = W_or_q0_dot';
%     W_or_q0_dot = NaN;   
% else
%     disp('Incorrect dimensions of either W_or_q0_dot or Ja')
% end

ra = [r_dot; y_dot];

if (nargin > 4)
    [mt, nt] = size (W_or_q0_dot);
    if (nt == 1)
        q0_dot = W_or_q0_dot;
        W_or_q0_dot = NaN;
        q_dot = LQO_sol(Ja, q0_dot, eye(n), ra);
    elseif (nt > 1)
        W = W_or_q0_dot;
        W_or_q0_dot = NaN;
        Jw_a = W_inv * jac';
        Jw_b = inv(jac * W_inv * jac');
        Jw = Jw_a * Jw_b;
        q_dot = Jw * ra;
    end
elseif (nargin == 4)
    q_dot = pinv(Ja) * ra;
else
    disp('Not enough input arguments');
end