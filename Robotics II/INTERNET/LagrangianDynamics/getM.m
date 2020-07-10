function M = getM (dh_table, sigma, q_dot, r_c, mass, I)
%USAGE: M = T(dh_table, sigma, qd, r_c, mass, I);
%dh_table in matrix of the form [alpha  a   d   theta]
%sigma = column vector of size = no of links; each element is zero if
%revolute and 1 otherwise
%q_dot = [qd1; qd2; qd3]
%r_c = [rc1, rc2, rc3]
%mass = [m1; m2; m3]; column vector
%I_c_zz(:,:,1) = I1, %I_c_zz(:,:,2) = I2, %I_c_zz(:,:,3) = I3
%I1 = diag([Ixx1, Iyy1, Izz1]);
%I2 = ...

% Robotics2Tools: https://github.com/MoSamArafat/Robotics2Tools
% Scripts based on Robotics 2 course at DIAG, Sapienza University of Rome taught by Prof. Alessandro De Luca.
% They are built to cover as many cases as found on previous exams.

% Most of the scripts have been tested on at least one example from previous exams.
% You are encouraged to test and provide feedback in case of any issues
 
% Contributors:
% Mohamed Samir Arafat
% Sayo Makinwa

[m,~] = size(dh_table);
w_prev = [0;0;0];
v_prev = [0;0;0];
KE = 0;

for i=1:m
    alpha = dh_table(i,1);
    a = dh_table(i,2);
    dM_dq = dh_table(i,3);

    theta = dh_table(i,4);
    rot_mat = [cos(theta)   ( -(cos(alpha)*sin(theta) ))    (sin(alpha)*sin(theta));
               sin(theta)   ( (cos(alpha)*cos(theta) ))     ( -(sin(alpha)*cos(theta)) );
               0            sin(alpha)                      cos(alpha)              ];

    r_i_i_1 = simplify([a*cos(theta); a*sin(theta); dM_dq]);

    %calculating T_i
    w_i_i_1 = simplify( w_prev + ( (1-sigma(i,1))*q_dot(i,1) ) * [0;0;1] );
    w_ii = simplify ( rot_mat' * w_i_i_1 )

    r_ii = simplify(rot_mat' * r_i_i_1);
    a = simplify(rot_mat' * ( v_prev + ( ( sigma(i,1)*q_dot(i,1) ) * [0;0;1] ) ));
    v_ii = simplify(a + cross(w_ii, r_ii) )

    vci = simplify ( v_ii + (cross(w_ii, r_c(:,i))) )

    T_i =  simplify ( ( 0.5 * mass(i,1) * (vci' * vci ) ) + (0.5 * w_ii' * I(:,:,i) * w_ii) );
    disp("T_i");
    pretty((T_i))
    disp("========================================================================");

    v_prev = v_ii;
    w_prev = w_ii;
    KE = KE + T_i;
end

KE = simplify(KE);

%Building the inertia matrix
for i=1:m
    q_d(i) = sym(q_dot(i));
end

%L = KE;
dLdq = sym(zeros(m,m));

for i = 1:m
    for k = 1:m
        dLdq(i,k) = diff(KE,q_d(1,i));
        dLdq(i,k) = simplify(dLdq(i,k));

        dLdq(i,k) = diff(dLdq(i,k),q_d(1,k));
        dLdq(i,k) = simplify(dLdq(i,k));
    end
end
M = simplify(dLdq);

% C terms
%     disp("C Terms");
%     disp("========================================================================");
%     C = sym(zeros(m,1));
%     S = sym(zeros(m,m));
%     for i=1:m
%         dM_dq = simplify(jacobian(M(:,i), q))
%         display(dM_dq')
%         display(diff(M,q(i)))
%         c_temp = simplify((0.5 .* (dM_dq + dM_dq' - diff(M,q(i)))))
%         S(:,i) = simplify( c_temp * q_dot);
%         C(i,1) = simplify( q_dot' * c_temp * q_dot )
%     end
%     
%     %To get S2, do C/qd, then make a comment that this solution is not
%     %unique. Then check if (ans * qd) is symmetric to verify if this usable
%     S2 = C/q_dot
%     
%     dM = sym(zeros(m,m));
%     for i = 1:m
%         for k = 1:m
%             for j = 1:m
%                 dM(i,k) = simplify (dM(i,k) + ((diff(M(i,k), q(j,1) )) * q_dot(j,1)));
%             end
%         end
%     end
%     dM
%     
%     symm1 = simplify ( dM - (2 * S) )
%     symm2 = simplify ( dM - (2 * S2) )
end