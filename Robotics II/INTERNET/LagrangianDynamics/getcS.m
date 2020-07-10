function [c, S] = getcS (M, q, q_dot)

% Robotics2Tools: https://github.com/MoSamArafat/Robotics2Tools
% Scripts based on Robotics 2 course at DIAG, Sapienza University of Rome taught by Prof. Alessandro De Luca.
% They are built to cover as many cases as found on previous exams.

% Most of the scripts have been tested on at least one example from previous exams.
% You are encouraged to test and provide feedback in case of any issues
 
% Contributors:
% Mohamed Samir Arafat
% Sayo Makinwa

[m, ~] = size(M);
c = sym(zeros(m,1));
S = sym(zeros(m,m));
sines = sym(zeros(1,m));
cosines = sym(zeros(1,m));

for i= 1:m
	sines(i) = sin(q(i));
    cosines(i) = cos(q(i));
end

vars = [sines ,cosines]

for i=1:m
    dM_dq   = collect(simplify(jacobian(M(:,i), q)), vars)
    c_i     = collect(simplify((0.5 .* (dM_dq + dM_dq' - diff(M,q(i))))), vars)
    c(i,1)  = collect(simplify(q_dot' * c_i * q_dot), vars);
    S(i,:)  = collect(simplify(c_i * q_dot), vars);
end
end