function [M, C, S] = Inertia(T,q,q_dot)
   
    [m,~] = size(q);
   
    
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
    
    %C terms
    disp("C Terms");
    disp("========================================================================");
    C = sym(zeros(m,1));
    S = sym(zeros(m,m));
    for i=1:m
        temp = simplify( jacobian(M(:,i), q) );
        c_temp = (0.5 * ( temp + temp' - diff(M,q(i)) ) );
        c_temp = simplify(c_temp);
        S(:,i) = simplify( c_temp * q_dot);
        C(i,1) = simplify( q_dot' * c_temp * q_dot );
    end
    
    %To get S2, do C/qd, then make a comment that this solution is not
    %unique. Then check if (ans * qd) is symmetric to verify if this usable
    S2 = C/q_dot
    
    dM = sym(zeros(m,m));
    for i = 1:m
        for k = 1:m
            for j = 1:m
                dM(i,k) = simplify ( dM(i,k) + diff(M(i,k), q(j,1)) );
            end
        end
    end
    dM
    
    symm1 = simplify ( dM - (2 * S) )
    symm2 = simplify ( dM - (2 * S2) )
end