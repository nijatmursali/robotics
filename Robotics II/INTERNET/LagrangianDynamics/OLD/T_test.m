syms L1 q q2 L2 real
dh_table = [-pi/2,  L1, 0,  q1;
             0      0   q2  0 ];
q = [ q1; q2];


sigma= [ 0 1 ]';

syms q1_d q2_d real
q_dot = [ q1_d; q2_d];

syms m1 m2 real
mass = [m1;m2];

syms dc1 dc2 real
rc_1 = [-L1+dc1; 0; 0];
rc_2=[ 0; 0 ; -dc2];

r_c=[ rc_1, rc_2]

syms I1 I2 I3 real
syms Ixx1 Ixx2 Iyy1 Iyy2 Izz1 Izz2 real
I1 = diag([Ixx1, Iyy1, Izz1]);
I2 =diag([Ixx2, Iyy2, Izz2]);
I_c_zz(:,:,1) = I1;
I_c_zz(:,:,2) = I2;


[M, C, S] = T(dh_table, sigma, q, q_dot, r_c, mass, I_c_zz)