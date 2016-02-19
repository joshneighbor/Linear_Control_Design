%%
fprintf('\nJosh Neighbor \nA53105202 \n280B');
close all; clear all; warning('off','all');

display('Problem 1');
display('part a');
s = tf('s');
display('we have our transfer function:');
G = (1.151*s + 0.1774) / (s^3 + .739*s^2 + .9215*s) % tf of plant G(s)

Poles = eig(G);% Check poles to see if stable

display('I looked at the root locus:')
figure; rlocus(G); title('Root Locus G(s)');
k = 3; %chose k based on root locus
fprintf(['I chose a gain of %.f which stabilized the system although this \n',...
    'would not meet specific design criteria for the dynamics of an aircraft',...
    '.\n'], k);
G_cl = G / ( 1 + G*k);
display('To check that my chosen gain stabilized the system I looked at the step responce of the closed loop system')
figure; step(G_cl); xlabel('time (sec)'); ylabel('pitch angle (rad)');
title('Closed-loop Step Response G(s)');
fprintf('part b \n Integrate G(s)')
G_hat = G / s
display(' Used root locus again to find a stabilizing value of K.')
figure; rlocus(G_hat)
k_hat = .3;
fprintf(['I chose a gain of %.f which stabilized the system.'], k_hat);
G_hat_cl = G_hat / ( 1 + G_hat * k_hat);
display('I then checked stability using step responce')
figure; step(G_hat_cl); xlabel('time (sec)'); ylabel('pitch angle (rad)');
title('Closed-loop Step Response $$\hat{G}$$ (s)','Interpreter','Latex');
display('part c')
fprintf(['The presence of the integrator in the plant guarantees zero steady-state error. \n']);
fprintf('Also introduces phase lag. An integrator is a form of first-order low-pass filter.\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
fprintf('\nProblem 2 LQR Design \nPart a \n')

A = [-.313 56.7 0; -.0139 -.426 0; 0 56.7 0];
B = [.232; .0203; 0]; C = [0 0 1]; D = [0];

fprintf(['I first made a state space system from A,B,C,D. I then took that\n'...
    ' system and make it into a transfer fucntion to check: G_sys_tf = tf(G_sys)\n']);
G_sys = ss(A, B, C, D); % creating system
G_sys_tf = tf(G_sys) % Matches our G(s) tf so this is the state space representation
display('We see that it matches the transfer fucntion in problem 1');
fprintf('\nPart b: controllable and observable?\n')
fprintf('I first built my controllability matrix CO = ctrb(A,B) and check the rank\n')
CO = ctrb(A,B);
rank_CO = rank(CO)
display('Note that I have full rank so is controllable. Now check rank of OB = obsv(A,C)');
OB = obsv(A,C);
rank_OB = rank(OB)
fprintf('note it is observable as well.\n\n')

display('Part c:  Design a feedback controller of form  u=Kx');
fprintf('\n\n\n\n\n\n\n');
fprintf('test for q = 1,10,100,1000\n')
q = [1 10 100 1000]; R = [1]; 
fprintf('Solved for K by using LQR. This solves the riccati equation by minimizing cost J.\n')
fprintf('\nPart d: \n')
for j = 1:4
    z = [0 q(j) 1; 1 0 0];
    Q_mat = z'*z;
    str = strcat('q', num2str(j));
    Q.(str) = Q_mat;
    [K_vec, S_vec, E_vec] = lqr(A, B, Q_mat, R);
    K.(str) = -K_vec; E.(str) = E_vec; S.(str) = S_vec;
end

fprintf('Calculating optimal cost when x0 = [1; 1; 1] by solving: \n')
fprintf('J = x0''*X*x0\n')
x0 = [1 1 1]';
for i = 1:4 
    num = strcat('q', num2str(i));
    J.(num) = x0'*S.(num)*x0;
end
fprintf('for a q value of 1, 10, 100, 1000 respectively we have a J of: \n')
J
fprintf('for a q of 1, 10, 100, 1000 respectively we have a K of: \n')
K
fprintf('for a q of 1, 10, 100, 1000 respectively we have closed loops poles (E): \n')
E.q1
E.q2
E.q3
E.q4


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
fprintf('\n\nProblem 3\n')
display('Part a & b');
fprintf('\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n');
display('Part c: Compute the open-loop transfer function H(s) for each of the LQR controllers');

for i = 1:4 
    num = strcat('q', num2str(i));
    fprintf('For a q value of %.f we have a H(s) of: \n', q(i))
    H.(num) = minreal( (-K.(num))*inv((s*eye(3)-A))*B);
    H.(num)
end

fprintf('\nPart D\n\nWe solve for each value of q: R(s) = H(s) / (1 + H(s))\n')

for i = 1:4 
    num = strcat('q', num2str(i));
    fprintf('For a q value of %.f we have a R(s)in ss of: \n', q(i))
    R.(num) = H.(num) / ( 1 + H.(num));
    R_sys.(num) = ss(R.(num))
   
end
fprintf('\nI also looked at the closed loop step responce for R(s) for q=1, attached\n')
figure; step(R.q1); ylabel('pitch angle (rad)');
title('Closed-Loop Step Response R(s) for q = 1: LQR');

fprintf('\nPart e\n');
fprintf(['Looking at R(s) = H(s) / (1 + H(s)) we can see that for a H(s) = -1\n'...
    'the system goes unstable and thus we can see the GM = 0DB \nand Phase = -180\n'...
    '\nPart f\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n']) % draw in bode example
fprintf(['Gain Margin:\n1)Find frequency where phase hits 180 degrees.\n'...
    '2) Look at Gain at this frequency.\n3) GM is 0 - gain@180.\nPhase Margin:'...
    '1) Find frequency where gain is 0 db.\n2) Find the phase at this freq.\n'...
    '3) PM = phase@0db + 180degrees.\n\nPart g\n\n'])
for i = 1:4 
    num = strcat('q', num2str(i));
    fprintf('For a q value %.f we have a H(s) w/ a GM of: \n', q(i))
    GM = margin(H.(num))
    figure; margin(H.(num))
end

fprintf('\nPart f\nGain Margin from Problem 1, where K = 3: \n')
GM_p1 = margin(G_hat_cl)
fprintf(['Note from our K chosen in problem 1 we do not have a good GM\n\n'...
    'Part i\nPhase Margin is the amount the phase would have to change for '...
    'the system to go unstable.\nThe further the phase is from 180 degress'...
    ' when the gain is 0bd the more robust our system is.\n This means that'...
    ' if our plant is not modeled perfectly we have some leeway.\n'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
fprintf('\nProblem 4 attached\n\nProblem 5    construct an estimator\n')
A5 = [-.313 56.7 0; -.0139 -.426 0; 0 56.7 0];
B5u = [.232; .0203; 0]; C5y = [0 0 1]; D5 = [0 0 1];
B5w = [.9 .4 0; .5 -.1 0; 1 -1.1 0];
eta = [1 10 100 1000];
fprintf('Part a\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n')%write in by hand


fprintf('\nPart b\n used lqr(A'',Cy'',Bw*W*Bw'',D*W*D'') to find my F for estimator.\n')
for i = 1:4
    var = strcat('eta', num2str(i));
    W.(var) = [eta(i) 0 0;0 eta(i) 0; 0 0 1]; %noise variance
    Q5 = B5w*W.(var)*B5w';
    R5 = D5*W.(var)*D5';
    [F_vec, Y_vec, E5_vec] = lqr(A5', C5y', Q5, R5);
    F.(var) = -F_vec'; E5.(var) = E5_vec'; Y.(var) = Y_vec';
    J5.(var) = trace(Y_vec);
    fprintf(['for the eta value of %.f we have an Optimal Cost (J), an\n'...
        'Optimal F matrix, and the poles of A+FC at the following: \n'], eta(i));
    J_mat = trace(Y_vec)
    F_mat = F_vec
    Poles = E5_vec
end
E_sys = reg(G_sys,-K.q3,-F.eta4); % to check with my OP_ctr

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
Cz = [0 100 1; 1 0 0];
display('Problem 6'); % optimal controller LQG
display('Part a   q = 100 for this question');
fprintf(['\nI took my K values from my solving the lqr, combined with my F'...
    ' values from my \nestimator to get my LQG controller.\n'...
    'OP_ctr = ss(A + Bu*K + F*Cy, -F, K, 0)\nAssumptions:\n1) (A,Bu) stabilizable\n'...
    '2) (A, Cy) is detectable\n3) w(t) is white with variance W>0\n'...
    '4) Cz''*Dz = 0 and Dz''*Dz > 0\n5) Bw*W*Dy''=0 and Dy*W*Dy''>0\n'])
fprintf('\nPart b\nI calculated cost using equation 6 from lecture 5 notes')
for j = 1:4
    var = strcat('eta', num2str(j));
    % building controller
    OP_ctr_sys = ss(A5 + B5u*K.q3 + F.(var) * C5y, F.(var), K.q3, 0);
    OP_ctr.(var) = OP_ctr_sys; % getting state space
    OP_ctr_tf_p = tf(OP_ctr.(var)); % getting tf
    OP_ctr_tf.(var) = OP_ctr_tf_p;
    fprintf(['For the eta value of %.f we have a sys and tf of: \n'], eta(j));
    OP_ctr_sys
    OP_ctr_tf_p
    H_new_new.(var) = ss([A5 zeros(3);F.(var)*C5y A5-F.(var)*C5y], [B5u; zeros(3,1)],[zeros(1,3) K.q3],0);
    H_new_new_tf.(var) = tf(H_new_new.(var));
    %Solve lyop
    Y1.(var) = lyap(A5-F.(var)*C5y, (B5w-F.(var)*D5)*W.(var)*(B5w-F.(var)*D5)');
    Y3.(var) = lyap(A5-B5u*K.q3, F.(var)*D5*W.(var)*D5'*F.(var)');
    %compute cost
    J6_vec = trace([Cz Cz]*[Y1.(var) zeros(3); zeros(3) Y3.(var)]*[Cz'; Cz']);
    J6.(var) = J6_vec;
    fprintf(['For the eta value of %.f we have a cost J of: \n'], eta(j));
    J6_vec
    % compute poles of closed loop system
    Acl.(var) = [A5, B5u * K.q3; -F.(var) * C5y, A5 + B5u * K.q3 + F.(var) * C5y];
    fprintf(['For the eta value of %.f we have closed loop poles at the following: \n'], eta(j));
    eig(Acl.(var))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
fprintf('\nProblem 7\nPart a\n')
display('Took the tf from our optimal controller in problem 6 and multiplied')
display('it by our G(s) to get our new open loop H(s). Noting that u_hat = -u.')
for j = 1:4
    var = strcat('eta', num2str(j));
    H_new_ss.(var) = ss(OP_ctr_tf.(var) * G);
    fprintf(['For the eta value of %.f we have our new H(s) in statespace: \n'], eta(j));
    H_new_ss.(var)
end
display('Part b')
for j = 1:4
    var = strcat('eta', num2str(j));
    H_new.(var) = OP_ctr_tf.(var) * G;
    fprintf(['For the eta value of %.f we have our new H(s): \n'], eta(j));
    H_new.(var)
end

display('Part c')
for i = 1:4 
    num = strcat('eta', num2str(i));
    fprintf('For eta value %.f we have our new H(s) w/ a GM of (not in bd): \n', eta(i))
    GM = margin(H_new.(num))
    figure; margin(H_new.(num))
end
fprintf(['\nWe see that with the introduction of noise our gain margins are'...
    ' no longer infinite.\n'])

display('Part d')
display('First I selected the noise varience with the best margin, which was an eta of 1000')
fprintf('\nCombined the new LQR gain with my estimator to form my new LQG controllers.\n')
sigma = [1e-1 1e-3 1e-6 1e-8 1e-10 1e-12 1e-14]; Q7 = C5y'*C5y;

for i = 1:7
    R = sigma(i)*eye(1);
    num = strcat('sigma', num2str(i));
    fprintf('for sigma value %1.0e we have our new LQR with K of: \n', sigma(i))
    [K_new_vec X_new_vec E_new_vec] = lqr(A5,B5u,Q7,R);
    K_new.(num) = -K_new_vec; E_new.(num) = E_new_vec; X_new.(num) = X_new_vec;
    -K_new_vec
end

for i = 1:7
   num = strcat('sigma', num2str(i));
   Est_sys.(num) = reg(G_sys,-K_new.(num),-F.eta4); %building lqg
end
fprintf('\nPart e\n')
for i = 1:7
   num = strcat('sigma', num2str(i));
   Est_sys_tf.(num) = tf(Est_sys.(num));
   fprintf('for sigma value %1.0e we have a Gain Margin (not in db) of: \n', sigma(i))
   GM = margin(-Est_sys.(num)*G)
   figure; margin(-Est_sys.(num)*G)
end
fprintf('\nPart f   Loop transfer recovery: gained / lost\n')
fprintf(['This method of finding tweaking our Q and R values we are able to'...
    ' get some of our \ngain margin back but at a greater cost.'])
fprintf('\nPart G\n')

%LQR design
%G_hat_sys = ss(G_hat); R7 = [1];
%for j = 1:4
    %str = strcat('q', num2str(j));
    %[K7_vec, S7_vec, E7_vec] = lqr(G_hat, Q7, R7);
    %K7.(str) = -K7_vec; E7.(str) = E7_vec; S7.(str) = S7_vec;
%end
%
