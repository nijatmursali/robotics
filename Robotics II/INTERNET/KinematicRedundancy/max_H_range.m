function H_range = max_H_range (q_curr, q_range)
%q_range is an m x 2 Matrix, where its first column is a vector of
%lower joint limits, and its second column is a vector of upper joint
%limits. q_range = [q_lower_limit, q_upper_limit];

% Robotics2Tools: https://github.com/MoSamArafat/Robotics2Tools
% Scripts based on Robotics 2 course at DIAG, Sapienza University of Rome taught by Prof. Alessandro De Luca.
% They are built to cover as many cases as found on previous exams.

% Most of the scripts have been tested on at least one example from previous exams.
% You are encouraged to test and provide feedback in case of any issues
 
% Contributors:
% Mohamed Samir Arafat
% Sayo Makinwa

[m, ~] = size (q_curr);
H_range = zeros(m,1);
q_bar = ((q_range(:,1) + q_range(:,2))./2);
q_temp = (q_range(:,2) - q_range(:,1));

SIGMA = 0;
for (i = 1:m)
    SIGMA = (((q_curr(i,1)-q_bar(i,1)) / (q_temp(i,1)))^2);
    H_range(i,1) = (1/(2*m)) * SIGMA;
end
H_range = double(H_range);

end