function M_dot = getM_dot (M, q, q_dot)

% Robotics2Tools: https://github.com/MoSamArafat/Robotics2Tools
% Scripts based on Robotics 2 course at DIAG, Sapienza University of Rome taught by Prof. Alessandro De Luca.
% They are built to cover as many cases as found on previous exams.

% Most of the scripts have been tested on at least one example from previous exams.
% You are encouraged to test and provide feedback in case of any issues
 
% Contributors:
% Mohamed Samir Arafat
% Sayo Makinwa

[m, ~] = size(M);
M_dot = sym(zeros(m,m));

for i = 1:m
    for j = 1:m
        for k = 1:m
            M_dot(i,j) = simplify (M_dot(i,j) + ((diff(M(i,j), q(k,1))) * q_dot(k,1)));
        end
    end
end

end