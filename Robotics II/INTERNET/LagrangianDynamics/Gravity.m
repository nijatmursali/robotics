function g_q = Gravity(dh_table, g_vector, mass, r_c, q)

% Robotics2Tools: https://github.com/MoSamArafat/Robotics2Tools
% Scripts based on Robotics 2 course at DIAG, Sapienza University of Rome taught by Prof. Alessandro De Luca.
% They are built to cover as many cases as found on previous exams.

% Most of the scripts have been tested on at least one example from previous exams.
% You are encouraged to test and provide feedback in case of any issues
 
% Contributors:
% Mohamed Samir Arafat
% Sayo Makinwa

clc
[m,~] = size(dh_table);
U = sym(zeros(m,1));
U_sum = 0;
T_mat = eye(4,4);

for i=1:m
    alpha = dh_table(i,1);
    a = dh_table(i,2);
    d = dh_table(i,3);

    theta = dh_table(i,4);

    T_mat_temp = [cos(theta)   ( -(cos(alpha)*sin(theta) ))    (sin(alpha)*sin(theta))      a*cos(theta);
               sin(theta)   ( (cos(alpha)*cos(theta) ))     ( -(sin(alpha)*cos(theta)) )    a*sin(theta);
               0            sin(alpha)                      cos(alpha)                      d;
               0            0                               0                               1 ];

    T_mat = simplify (T_mat * T_mat_temp);

    rot_mat = T_mat(1:3,1:3);               %Concatinated Rotations Up to Link i
    ri = T_mat(1:3,4);                      %Distance from Base Frame to Origin of Link i expressed in the Base frame
    rc_0i = simplify(rot_mat * r_c(:,i));   %Distance from the Origin of Link i to the Center of Mass of Link i
    rc = simplify (ri + (rc_0i));           %Distance from the Base Frame to the Center of Mass of Link i, expressed in the Base Frame
    U(i,1) = simplify (- mass(i,1) * g_vector' * rc )
    U_sum = simplify (U_sum + U(i,1));
end

g_q = jacobian(U_sum, q)';
end