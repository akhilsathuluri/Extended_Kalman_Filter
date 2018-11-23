tic
%This code takes in data and gives better estimates using Kalman Filter
load('DATA_REQD_2.mat')
%Variables
ts = 17.6039;
I = eye(3);
%Intiating the state and co-variance matrices
x0 = [0.00385;0.00472;393]; %Assuming the initial values of the other state to be thier nominal values
P0 = [10 0 0;0 10 0;0 0 100];
% P0 = 10*rand(3);
%The error covariances of the system model and the measurement model are
%assumed here
Q = [1 0 0;0 1 0;0 0 1]*10;
% R = [1 0 0;0 1 0;0 0 1];
R = 10^3;
%Since only one output is measured
Cmat = [0 0 1];
xhk_k = x0;
Pk_k = P0;
i = 1;
% J = JacobianSystem;

for i = 1:(size(Y_measured_case1)-1)
% tempJ = J;
%% Prediction step of the states, by solving the ode model for a given time step
[tp,xp] = ode45(@(t,x)fcc_fn_to_solve_odemodel(t,x,InputU(i,:)),[Time(i),Time(i+1)],[xhk_k]);
xhk1_k = transpose(xp(end,:));
%% Prediction step of the covariance matrix
%For predicting the covariance matrix we'll need to linearize the system
%about the current state values.
% C_rc = xhk_k(1);   O_d = xhk_k(2);  T_rg = xhk_k(3); 
% dFa  = InputU(i,1); dFsc  = InputU(i,2);   
% tempA = double(subs(tempJ));
Amat = JacobianSystem(xhk_k, InputU(i,:));
% Amat = exp(tempA*ts);
Pk1_k = Amat*Pk_k*transpose(Amat)+Q;
%% Now computation of the Kalman gain
K = (Pk1_k*transpose(Cmat))/(Cmat*Pk1_k*transpose(Cmat)+R);
%% Correction step of state
xhk1_k1 = xhk1_k + K* (Y_measured_case1(i+1)-Cmat*xhk1_k);
%% Correction step of covariance matrix
Pk1_k1 = (I-K*Cmat)*Pk1_k;
%% Now updating all the older variables of the loop
xhk_k = xhk1_k1;
Pk_k = Pk1_k1;
xg(i,:) = xhk1_k1;
i = i+1;
end
plot(Time, Y_measured_case1, 'g');
hold on;
plot(Time, xg(1:2046,3),'r');  
title('Kalman Filter in Action');
xlabel('Time(s)');
ylabel('Measured Data');
legend({'Measured Data','Filtered Data'},'Location','southeast');

%plot(xg(:,3),'r');  
xg = real([transpose(x0);xg]);
toc