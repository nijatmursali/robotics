function K = UElastic(dh_table, g_vector, mass, r_c, q, theta, n, K)

% Robotics2Tools: https://github.com/MoSamArafat/Robotics2Tools
% Scripts based on Robotics 2 course at DIAG, Sapienza University of Rome taught by Prof. Alessandro De Luca.
% They are built to cover as many cases as found on previous exams.

% Most of the scripts have been tested on at least one example from previous exams.
% You are encouraged to test and provide feedback in case of any issues
 
% Contributors:
% Mohamed Samir Arafat
% Sayo Makinwa

clc
[m,~] = size(q);

Ue = 0;
for i=1:m
    Uei = simplify( 0.5 * K(i,i) * ( ( q(i,1) - ( theta(i,1)/n(i,1) ) )^2 ) )
    Ue = simplify( Ue + Uei );
end

g_q = Gravity(dh_table, g_vector, mass, r_c, q)

K = sym(zeros(m,2));

K(:,1) = jacobian(Ue, q)';
K(:,2) = jacobian(Ue, theta)';
end