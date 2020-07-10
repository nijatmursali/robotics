function q_dot = RG_sol(jac, q0_dot, r_dot)
% 02_KinematicRedundancy.pdf, page 33*
% conversion to rad from deg: deg*pi/180
% jac = double(jac);

% Robotics2Tools: https://github.com/MoSamArafat/Robotics2Tools
% Scripts based on Robotics 2 course at DIAG, Sapienza University of Rome taught by Prof. Alessandro De Luca.
% They are built to cover as many cases as found on previous exams.

% Most of the scripts have been tested on at least one example from previous exams.
% You are encouraged to test and provide feedback in case of any issues
 
% Contributors:
% Mohamed Samir Arafat
% Sayo Makinwa

numbers = [1, 2, 3, 4, 5, 6, 7, 8, 9];
[m, n] = size(jac);
numbers = numbers(1, 1:n);
ordera = zeros(1,m);
orderb = zeros(1,n-m);
rho = rank (jac);
DET = 0;
maxDET = 0;

%Finding the non-singular minor with the largest determinant:
if (rho == m && m < n)
    combos = nchoosek(1:n,m);
    [mc,~] = size(combos);
    minor_jac = zeros(m, m);
    Ja = zeros (m,m);
    Jb = zeros(m, n-m); %(m x (n-m))

    for (i = 1:mc)
        for j = 1:m
            minor_jac(:,j) = jac(:,combos(i,j));
            ordera(j) = numbers(1, combos(i,j));
        end
        DET = abs(det(minor_jac));
        if (DET >= maxDET)
            maxDET = DET;
            if (DET > 0.01)
                Ja = minor_jac;
            end
        end
    end
    
    Jb = setdiff(jac',Ja','rows');
    Jb = Jb';
    display(Ja)
    display(maxDET)
    display(Jb)
    
    [ma, na] = size (Ja);
    [mb, nb] = size (Jb);
    [mq0, nq0] = size (q0_dot);
    
    %Calculating q_dot!
    Ja_inv = inv(Ja);
    Q1 = [Ja_inv; zeros(mq0-ma,na)] * r_dot;
    Q2 = [-1*(Ja_inv*Jb); eye(mq0-ma,nb)];
    Q3 = [-1*((Ja_inv*Jb)'), eye(nb, mq0-ma)];
    q_dot_t = Q1 + (Q2*Q3*q0_dot);
    
    orderb = setdiff(numbers, ordera);
    order = fliplr([ordera, orderb]);
    for (i = 1:n)
    	q_dot(i,1) = q_dot_t(order(1,i));
    end
else
    disp("Rank and Size conditions are not met!")
    return
end
end