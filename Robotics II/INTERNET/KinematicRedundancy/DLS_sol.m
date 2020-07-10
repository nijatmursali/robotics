function [q_dot, q_dot_norm, err_dot_norm] = DLS_sol(jac, r_dot, mu_squared)
%function [q_dot, q_dot_norm, err_dot_norm, e_dot, H] = DLS_solution(jac, r_dot, mu_squared)
%02_KinematicRedundancy.pdf, page 24* - 26
%conversion to rad from deg: deg*pi/180

% Robotics2Tools: https://github.com/MoSamArafat/Robotics2Tools
% Scripts based on Robotics 2 course at DIAG, Sapienza University of Rome taught by Prof. Alessandro De Luca.
% They are built to cover as many cases as found on previous exams.

% Most of the scripts have been tested on at least one example from previous exams.
% You are encouraged to test and provide feedback in case of any issues
 
% Contributors:
% Mohamed Samir Arafat
% Sayo Makinwa

jac = double(jac);
[m, n] = size(jac);
[m_, ~] = size(r_dot);

if (m == m_)
    [U, S, V] = svd(jac);
    rho = rank(jac);

    S_dls1 = zeros(rho, rho);
    S_dls2 = zeros(m-rho, m-rho);
    S_dls3 = zeros(n-m, rho);
    S_dls4 = zeros(n-m, m-rho);

    for (i = 1:rho) 
        S_dls1(i,i) = S(i,i) / ((S(i,i)^2) +  mu_squared);
    end
    S_dls = [S_dls1 S_dls2; S_dls3 S_dls4];
    J_dls = V*S_dls*U';
    q_dot = J_dls*r_dot;
    q_dot_norm = norm(q_dot);
    e_dot = mu_squared * (inv(jac*jac' + mu_squared * eye(m))) * r_dot;
    H = ((mu_squared/2) * (q_dot_norm^2)) + (0.5*((norm(r_dot - jac*q_dot))^2));
    err_dot_norm = norm(e_dot);
else
    display("You need to clip the Jacobian or something")
    return
end
end