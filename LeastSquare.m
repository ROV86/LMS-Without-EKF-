% Initialization
time_steps = 100;         % Number of time steps in the simulation
dt = 1;                   % Time step size
R = (0.01)^2;             % 1cm/s DVL standard deviation

M = [25, 25, 25, 10, 10, 10]'; % Mass-Inertia Matrix

% Initial Estimates
D_real = [0.35; 0.5; 0.35; 1; 0.4; 0.4]; % Improved initial guesses for drag coefficients
mA_real = [15; 15; 10; 5; 8; 5]; % Real added masses


% Thruster Allocation Matrix
T = [0.707  0.707   -0.707  -0.707  0   0   0   0;
     -0.707 0.707   -0.707  0.707   0   0   0   0;
     0  0   0   0   -1  1   1   -1;
     0.06   -0.06   0.06   -0.06  -0.218  -0.218  0.218  0.218;
     0.06   0.06   -0.06  -0.06   0.120  -0.120  0.120  -0.120;
    -0.188 0.188   0.188 -0.188  0   0   0   0]; % Defines how thruster forces contribute to system's degrees of freedom

% Constant Control Input
tau_input = [5; 5; -5; -5; 0; 0; 0; 0]; % Control input applied to the system

% Simulated dynamics and measurements
v_measured = zeros(6, time_steps); % Measured velocities over time

% Simulation loop
v_measured(:,1) = 0;

for k = 2:time_steps
    % Calculate thruster forces
    tau_input = 5*randn(8,1);% Control input applied to the system
    tau = T * tau_input;       % Compute forces applied by the thrusters
    % Dynamics equations
    %a_p(:,k) = tau ./ (M + mA_k) - D_k .* v_measured(:, k-1)/(M+mA_k); % Predicted acceleration
    a_m(:,k) = tau ./ (M + mA_real) - D_real .* v_measured(:, k-1)./(M + mA_real);% real acceleration
    
    v_measured(:, k) = v_measured(:, k-1) + a_m(:,k) * dt + randn(6,1)*sqrt(R);; % Update real velocity
    %v_predicted(:, k) = v_measured(:, k-1) + a_p(:,k) * dt; % Update predicted velocity
   
    % Acceleration_measured
    a_m = (v_measured(:,k)-v_measured(:,k-1))/dt; 
    %Least Square equation for mutliple measurements
    if (k == 2)
        A = [diag(a_m) diag(v_measured(:,k-1))];
        b = [tau];
    else
        A = [A;[diag(a_m) diag(v_measured(:,k))]];
        b = [b;tau];
    end
    
end
X = lsqminnorm(A,b);
mA_k = X(1:6) - M
D_k = X(7:12)


