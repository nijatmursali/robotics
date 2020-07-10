function [f, t, up] = NE_method (dh_table, sigma, q_d, q_dd , r_ci,  g, mass, I)

% Robotics2Tools: https://github.com/MoSamArafat/Robotics2Tools
% Scripts based on Robotics 2 course at DIAG, Sapienza University of Rome taught by Prof. Alessandro De Luca.
% They are built to cover as many cases as found on previous exams.

% Most of the scripts have been tested on at least one example from previous exams.
% You are encouraged to test and provide feedback in case of any issues
 
% Contributors:
% Mohamed Samir Arafat
% Sayo Makinwa

[m, ~]  = size (dh_table);
syms g0 q real
q_sum = 0;
q = sym(zeros(1,m));
% w   = sym(zeros(3,m));
% w_d = sym(zeros(3,m));
% a   = sym(zeros(3,m));
% ac  = sym(zeros(3,m));

w_init      = [0 0 0]';
w_d_init    = [0 0 0]';
a_init      = -g;

%Initial Step in the Forward Pass

alpha   = dh_table(1,1);
adh     = dh_table(1,2);
d       = dh_table(1,3);
theta   = dh_table(1,4);

rot_mat = [cos(theta)   ( -(cos(alpha)*sin(theta) ))    (sin(alpha)*sin(theta));
           sin(theta)   ( (cos(alpha)*cos(theta) ))     ( -(sin(alpha)*cos(theta)) )
           0            sin(alpha)                      cos(alpha)              ];

r_i_i_1 = simplify([adh*cos(theta); adh*sin(theta); d]);

if (sigma(1,1) == 0)        %Revolute Joint
    disp("Joint 1 Revolute")
    q_sum = theta;
    q(1) = q_sum;
    w(:, 1)     = simplify(rot_mat' * (w_init + (q_d(1) * [0 0 1]')))
    w_d(:, 1)   = simplify(rot_mat' * (w_d_init + (q_dd(1)*[0 0 1]') + cross((q_d(1)*w_init), [0 0 1]')))
    a(:, 1)     = simplify((rot_mat' * a_init) + cross(w_d(:,1), (rot_mat' * r_i_i_1)) + cross(w(:,1), cross(w(:,1), (rot_mat' * r_i_i_1))));
%     disp("r_01")
%     disp(simplify(rot_mat' * r_i_i_1))
    a(:,1)      = collect(a(:,1), [q_dd(1) q_d(1) theta])
    ac(:, 1)    = simplify(a(:,1) + cross(w_d(:,1), r_ci(:,1)) + cross(w(:,1), cross(w(:,1), r_ci(:,1))));
    ac(:, 1)    = collect(ac(:, 1), [q_dd(1) q_d(1) theta])
elseif (sigma(1,1) == 1)    %Prismatic Joint
    disp("Joint 1 Prismatic")
    w(:, 1)     = simplify(rot_mat' * w_init)
    w_d(:, 1)   = simplify(rot_mat' * w_d_init)
    a(:, 1)     = simplify(rot_mat' * (a_init+(q_dd(1)*[0 0 1]')) + cross((2*q_d(1)*w(:,1)), (rot_mat'*[0 0 1]')) + cross(w_d(:,1), (rot_mat' * r_i_i_1)) + cross(w_init , cross(w_init, rot_mat' * r_i_i_1)))
    a(:,1)      = collect(a(:,1), [q_dd(1) q_d(1) theta])
%     disp("r_01")
%     disp(simplify(rot_mat' * r_i_i_1))
    ac(:, 1)    = simplify(a(:,1) + cross(w_d(:,1), r_ci(:,1)) + cross(w(:,1), cross(w(:,1), r_ci(:,1))));
    ac(:, 1)    = collect(ac(:, 1), [theta q_d(1) q_dd(1)])
end
disp("------------------------------------------------------------------------------------------------------------------")
%=========================================================================================================================

for (i=2:m)  %Forward Pass
    
    alpha   = dh_table(i,1);
    adh     = dh_table(i,2);
    d       = dh_table(i,3);
    theta   = dh_table(i,4);
    
    rot_mat = [cos(theta)   ( -(cos(alpha)*sin(theta) ))    (sin(alpha)*sin(theta));
               sin(theta)   ( (cos(alpha)*cos(theta) ))     ( -(sin(alpha)*cos(theta)) )
               0            sin(alpha)                      cos(alpha)              ];
    
    r_i_i_1 = simplify([adh*cos(theta); adh*sin(theta); d]);
    
    if (sigma(i,1) == 0)        %Revolute Joint
        disp("Revolute")
        q_sum = q(i-1) + theta;
        q(i) = q_sum;
        w(:, i)     = simplify(rot_mat' * (w(:,i-1) + (q_d(i) * [0 0 1]')))
        w_d(:, i)   = simplify(rot_mat' * (w_d(:,i-1) + (q_dd(i)*[0 0 1]') + cross((q_d(i)*w(:,i-1)), [0 0 1]')))
        a(:, i)     = simplify((rot_mat' * a(:,i-1)) + cross(w_d(:,i), (rot_mat' * r_i_i_1)) + cross(w(:,i), cross(w(:,i), (rot_mat' * r_i_i_1))));
        a(:,i)      = collect(a(:,i), [q_dd(i) q_d(i) theta g0])
%         disp("r_0i")
%         disp(simplify(rot_mat' * r_i_i_1))
        ac(:, i)    = simplify(a(:,i) + cross(w_d(:,i), r_ci(:,i)) + cross(w(:,i), cross(w(:,i), r_ci(:,i))));
        ac(:,i)     = collect(ac(:,i), [q_dd(i) q_d(i) theta g0])
    elseif (sigma(i,1) == 1)    %Prismatic Joint
        disp("Prismatic")
        w(:, 1)     = simplify(rot_mat' * w(:,i-1))
        w_d(:, 1)   = simplify(rot_mat' * w_d(:, i-1))
        a(:, 1)     = simplify(rot_mat' * (a(:,i-1)+(q_dd(i)*[0 0 1]')) + cross((2*q_d(i)*w(:,i)), (rot_mat'*[0 0 1]')) + cross(w_d(:,i), (rot_mat' * r_i_i_1)) + cross(w(:,i-1) , cross(w(:,i-1), rot_mat' * r_i_i_1)));
        a(:,i)      = collect(a(:,i), [q_dd(i) q_d(i) theta g0])
%         disp("r_0i")
%         disp(simplify(rot_mat' * r_i_i_1))
        ac(:, i)    = simplify(a(:,i) + cross(w_d(:,i), r_ci(:,i)) + cross(w(:,i), cross(w(:,i), r_ci(:,i))));
        ac(:,i)     = collect(ac(:,i), [q_dd(i) q_d(i) theta g0])
    end
disp("-------------------------------------  !END OF FORWARD RECURSION!  --------------------------------------------------------------")
end
%End of Forward Path
%===================================================================================================================

f_init = [0 0 0]';
t_init = [0 0 0]';
% f   = sym(zeros(3,m));
% t   = sym(zeros(3,m));
up  = sym(zeros(m,1));

%Initial Step in the Backward Path

alpha   = dh_table(m,1);
adh     = dh_table(m,2);
d       = dh_table(m,3);
theta   = dh_table(m,4);
rot_mat = [cos(theta)   ( -(cos(alpha)*sin(theta) ))    (sin(alpha)*sin(theta));
           sin(theta)   ( (cos(alpha)*cos(theta) ))     ( -(sin(alpha)*cos(theta)) )
           0            sin(alpha)                      cos(alpha)              ];

r_i_i_1 = simplify([adh*cos(theta); adh*sin(theta); d]);

if (sigma(m,1) == 0)        %Revolute Joint
    disp("Joint m: Revolute")
    f(:,m)  = simplify(f_init + (mass(m)* ac(:,m)));
    f(:,m)  = collect(f(:,m), [q_dd(m) q_d(m) theta g0])
    t(:,m)  = simplify(t_init - cross(f(:,m), ((rot_mat'*r_i_i_1)+r_ci(:,m))) + cross(f_init, r_ci(:,m)) + (I(:,:,m)*w_d(:,m)) + cross(w(:,m), (I(:,:,m)*w(:,m))));
    t(:,m)  = collect(t(:,m), [q_dd(m) q_d(m) theta g0])
    up(m)   = simplify(t(:,m)' * rot_mat' * [0 0 1]');
    up(m)   = collect(up(m), [q_dd(m) q_d(m) theta g0])
elseif (sigma(m,1) == 1)    %Prismatic Joint
    disp("Joint m: Prismatic")
    f(:,m)  = simplify(f_init + (mass(m)* ac(:,m)));
    f(:,m)  = collect(f(:,m), [q_dd(m) q_d(m) g0])
    t(:,m)  = simplify(t_init - cross(f(:,m), ((rot_mat'*r_i_i_1)+r_ci(:,m))) + cross(f_init, r_ci(:,m)) + (I(:,:,m)*w_d(:,m)) + cross(w(:,m), (I(:,:,m)*w(:,m))));
    t(:,m)  = collect(t(:,m), [q_dd(m) q_d(m) g0])
    up(m)   = simplify(f(:,m)' * rot_mat' * [0 0 1]');
    up(m)   = collect(up(m), [q_dd(m) q_d(m) g0])
end
disp("----------------------------------------------------------------------------------------------------------------------")    
for (i = (m-1):-1:1) %Backward path

    alpha   = dh_table(i,1);
    adh     = dh_table(i,2);
    d       = dh_table(i,3);
    theta   = dh_table(i,4);
    rot_mat_next = [cos(theta)   ( -(cos(alpha)*sin(theta) ))    (sin(alpha)*sin(theta));
                   sin(theta)   ( (cos(alpha)*cos(theta) ))     ( -(sin(alpha)*cos(theta)) )
                   0            sin(alpha)                      cos(alpha)              ];
               
    alpha   = dh_table(i+1,1);
    adh     = dh_table(i+1,2);
    d       = dh_table(i+1,3);
    theta   = dh_table(i+1,4);
    rot_mat = [cos(theta)   ( -(cos(alpha)*sin(theta) ))    (sin(alpha)*sin(theta));
               sin(theta)   ( (cos(alpha)*cos(theta) ))     ( -(sin(alpha)*cos(theta)) )
               0            sin(alpha)                      cos(alpha)              ];

    r_i_i_1 = simplify([adh*cos(theta); adh*sin(theta); d]);
    
    
    if (sigma(i,1) == 0)        %Revolute Joint
        disp("Revolute")
        f(:,i)  = simplify((rot_mat * f(:,i+1)) + (mass(i)* (ac(:,i))));
        f(:,i)  = collect(f(:,i), [q_dd' q_d sin(q) cos(q) g0])
        t(:,i)  = simplify((rot_mat * t(:,i+1)) - cross(f(:,i), ((rot_mat'*r_i_i_1) + r_ci(:,i))) + cross((rot_mat * f(:,i+1)), r_ci(:,i)) + (I(:,:,i)*w_d(:,i)) + cross(w(:,i), (I(:,:,i)*w(:,i))));
        t(:,i)  = collect(t(:,i), [q_dd' q_d sin(q) cos(q) g0])
        up(i)   = simplify(t(:,i)' * rot_mat_next' * [0 0 1]');
        up(i)   = collect(up(i), [q_dd' q_d sin(q) cos(q) g0])
        elseif (sigma(i,1) == 1)    %Prismatic Joint
        disp("Prismatic")
        f(:,i)  = simplify((rot_mat * f(:,i+1)) + (mass(i)* ac(:,i)));
        f(:,i)  = collect(f(:,i), [q_dd' q_d sin(q) cos(q) g0])
        t(:,i)  = simplify((rot_mat * t(:,i+1)) - cross(f(:,i), ((rot_mat'*r_i_i_1) + r_ci(:,i))) + cross((rot_mat * f(:,i+1)), r_ci(:,i)) + (I(:,:,i)*w_d(:,i)) + cross(w(:,i), (I(:,:,i)*w(:,i))));
        t(:,i)  = collect(t(:,i), [q_dd' q_d sin(q) cos(q) g0])
        up(i)   = simplify(f(:,i)' * rot_mat_next' * [0 0 1]');
        up(i)   = collect(up(i), [q_dd' q_d sin(q) cos(q) g0])
    end
end
disp("------------------------------------------------------------------------------------------------------------------")
end

