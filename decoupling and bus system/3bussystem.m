% Define system parameters
V = [1.06; 1.0; 1.0]; % Initialize voltage magnitudes
P = [0; 0.65; -3.42]; % Real power injections
Q = [0; 0; -0.96];    % Reactive power injections

% Admittance matrix (Ybus) for a 3-bus system
Y = [0+0j, -1/(0.018j), -1/(0.021j);
     -1/(0.018j), 0+0j, -1/(0.028j);
     -1/(0.021j), -1/(0.028j), 0+0j];

% Set convergence tolerance and maximum iterations
epsilon = 1e-6;
max_iterations = 50;
iteration = 0;
delta_V = ones(3, 1);

% Store convergence values for analysis
convergence = [];

while max(abs(delta_V)) > epsilon && iteration < max_iterations
    V_prev = V;
    J = zeros(3);

    for i = 2:3 % Start from bus 2 to avoid reference bus
        delta_P = 0;
        delta_Q = 0;

        for j = 1:3
            delta_P = delta_P + V(i) * V(j) * (real(Y(i, j)) * cos(angle(V(i)) - angle(V(j))) + imag(Y(i, j)) * sin(angle(V(i)) - angle(V(j))));
            if i ~= j
                delta_Q = delta_Q + V(i) * V(j) * (real(Y(i, j)) * sin(angle(V(i)) - angle(V(j))) - imag(Y(i, j)) * cos(angle(V(i)) - angle(V(j))));
            end
        end

        % Jacobian matrix calculation
        if i == 2
            J(i, i) = -2 * V(i) * (P(i) - delta_P) - 2 * V(i)^2 * real(Y(i, i));
            J(i, 1) = V(i) * (P(i) - delta_P) * sin(angle(V(i)) - angle(V(1))) - V(i)^2 * real(Y(i, 1)) * sin(angle(V(i)) - angle(V(1)));
            J(i, i+1) = -2 * V(i) * (P(i) - delta_P) - 2 * V(i)^2 * real(Y(i, i+1));
            J(i, i+2) = -2 * V(i) * (P(i) - delta_P) - 2 * V(i)^2 * real(Y(i, i+2));
        elseif i == 3
            J(i, i) = -2 * V(i) * (P(i) - delta_P) - 2 * V(i)^2 * real(Y(i, i));
            J(i, 1) = V(i) * (P(i) - delta_P) * sin(angle(V(i)) - angle(V(1))) - V(i)^2 * real(Y(i, 1)) * sin(angle(V(i)) - angle(V(1)));
            J(i, i-1) = -2 * V(i) * (P(i) - delta_P) - 2 * V(i)^2 * real(Y(i, i-1));
            J(i, i) = -2 * V(i) * (P(i) - delta_P) - 2 * V(i)^2 * real(Y(i, i));
        end
    end

    % Solve for the change in voltages
    delta_V = J \ [P(2) - delta_P; P(3) - delta_P; Q(3) - delta_Q];
    V(2:3) = V(2:3) + delta_V(1:2);
    iteration = iteration + 1;
    convergence = [convergence, max(abs(delta_V))];
end

% Display and analyze results
if iteration < max_iterations
    disp('Converged to a solution!');
    disp('Voltage magnitudes:');
    disp(abs(V));
    disp('Voltage angles:');
    disp(angle(V));
else
    disp('Did not converge within the specified iterations.');
end

% Display convergence characteristics
disp('Convergence Characteristics:');
disp(convergence);