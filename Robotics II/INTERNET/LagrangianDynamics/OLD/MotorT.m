%function [Tm] = MotorT(dh_table, sigma, q, q_dot, r_c, mass, I_c_zz, Im, theta, n)
function [Tm] = MotorT(dh_table, q, q_dot, mass, Im, n)
    %Page 13 5a_slides
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
    % theta_n could be n or theta according to the data
    %m_option= theta if we insert theta or other stuff if is n
   
    [m,~] = size(dh_table);
    
    T_m = 0;
    T_mi = 0;
    w_m_init = [0 0 0]'
    w_m = sym(zeros(3,m))
    p_d = sym(zeros(3,m))
    
    alpha = dh_table(1,1);
    a = dh_table(1,2);
    d = dh_table(1,3);
    theta = dh_table(1,4);
    T_mat = [cos(theta)   ( -(cos(alpha)*sin(theta) ))    (sin(alpha)*sin(theta))         a*cos(theta);
               sin(theta)   ( (cos(alpha)*cos(theta) ))     ( -(sin(alpha)*cos(theta)) )    a*sin(theta);
               0            sin(alpha)                      cos(alpha)                      d;
               0            0                               0                               1 ];
    rot_mat = T_mat(1:3,1:3)
           
    w_m(:,1)    = w_m_init + (n(1,:) * q_dot(1,:) * (rot_mat * [0 0 1]'))
    T_m1        = 0.5 * w_m(:,1)' * Im(:,:,1) * w_m(:,1)
    
    
    for (i= 2:1:m)
        
        alpha = dh_table(i,1);
        a = dh_table(i,2);
        d = dh_table(i,3);
        theta = dh_table(i,4);
        T_mat = T_mat * [cos(theta)  ( -(cos(alpha)*sin(theta) ))    (sin(alpha)*sin(theta))        a*cos(theta);
                            sin(theta)   ( (cos(alpha)*cos(theta) ))     ( -(sin(alpha)*cos(theta)) )    a*sin(theta);
                            0            sin(alpha)                      cos(alpha)                      d;
                            0            0                               0                               1 ];
        rot_mat = T_mat(1:3,1:3)
                        
        w_m(:,1)    = w_m(:,i-1) + (n(i,:) * q_dot(i,:) * (rot_mat * [0 0 1]'))
        p(:,i)      = T_mat (1:3,4)
        for (j = 1:m)
            p_d(:,i) = p_d(:,i) + jacobian(p(:,i),(q(j))) .* q_dot(j);
        end
        T_mi = 0.5 * (mass(i,1) * (p_d(:,i)' * p_d(:,i))) + 0.5 * (w_m(:,i)' * Im(:,:,i) * w_m(:,i))
        T_m = T_m + T_mi;  
    end
    Tm = T_m1 + T_m
    
    %Building the BM matrix
    
%     for i=1:m
%         q_d(i) = sym(q_dot(i));
%     end
%     
%     dBdq = sym(zeros(m,m));
%     
%     for i = 1:m
%         for k = 1:m
%             dBdq(i,k) = diff(T_m,q_d(1,i));
%             dBdq(i,k) = simplify(dBdq(i,k));
%             
%             dBdq(i,k) = diff(dBdq(i,k),q_d(1,k));
%             dBdq(i,k) = simplify(dBdq(i,k));
%         end
%     end
%     Bm = simplify(dBdq);
%     
%     [M, C, S] = T(dh_table, sigma, q, q_dot, r_c, mass, I_c_zz);
%     
end