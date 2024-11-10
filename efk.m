close all;
clear all;

%% Initializing 

no_of_samples = 100;
dt = 0.1;
R1 = 0.1; % Measurement noise covariance for x
R2 = 0.3; % Measurement noise covariance for y
R3 = 0.4; % Measurement noise covariance for theta
t = 0 : dt : no_of_samples * dt; 
V = 2;
omega = 0.5;
x_true = zeros(1, length(t));
y_true = zeros(1, length(t));
theta_true = zeros(1, length(t));
x_true(1) = 10;
y_true(1) = 20;
theta_true(1) = 10;
true_val = zeros(3, length(t));
true_val(:,1) = [x_true(1);
                 y_true(1);
                 theta_true(1)];
noise_val = zeros(3, length(t));
for i = 2:length(t)
    S = V * dt;
    x_true(i) = x_true(i-1) -S * cos(theta_true(i-1) + (omega * dt / 2));
    y_true(i) = y_true(i-1) + S * sin(theta_true(i-1) + (omega * dt / 2));
    theta_true(i) = theta_true(i-1) + omega * dt;
    true_val(:,i) = [x_true(i);
                     y_true(i);
                     theta_true(i)];
end
%% Adding Noise 

for i = 1:length(t)
    noise_val(:,i) = [x_true(i) + R1 * randn();
                      y_true(i) + R2 * randn();
                      theta_true(i) + R3 * randn()];
end
%% EKF Initialization

X = zeros(3, length(t));
X(:,1) = [5; 15; 0];  % Adjusted initial state estimate   
Q = diag([0.0019, 0.007764, 0.164]);
H = [1 0 0;
     0 1 0;
     0 0 1];  
P = [8.8^2 0    0 ; 
     0 5.4^2 0 ;
     0 0    20.9^2];
R = diag([R1^2, R2^2, R3^2]);  

for i = 2:length(t)

    S = V * dt;
    X(:,i) = [X(1,i-1) + -S*cos(X(3,i-1) + (omega*dt/2));
              X(2,i-1) + S*sin(X(3,i-1) + (omega*dt/2));
              X(3,i-1) + omega*dt];
    
  
    F = [1 0 -S*sin(X(3,i-1) + (omega*dt/2));
         0 1  S*cos(X(3,i-1) + (omega*dt/2));
         0 0  1];
    
 
    P = F*P*F' + Q;  
    

    K = P * H' / (H * P * H' + R);
    
 
    X(:,i) = X(:,i) + K * (noise_val(:,i) - H * X(:,i));
    

    P = (eye(3) - K * H) * P;
end

%% Plotting Results as a Simulation

for i = 1:length(t)

    clf;
    

   plot(x_true(1:i), y_true(1:i), 'b', 'LineWidth', 2); 
   hold on;
    

    plot(noise_val(1,1:i), noise_val(2,1:i), '.r');
    
    hold on;
    plot(X(1,1:i), X(2,1:i), 'g', 'LineWidth', 2);
    

    legend( 'ground truth','Noisy Measurements', 'Estimated Path');
    xlabel('X Position');
    ylabel('Y Position');
    title('EKF Estimation of Robot Trajectory');
    grid on;
    

    pause(0.01); 
    

    axis([min(x_true) max(x_true) min(y_true) max(y_true)]);
end