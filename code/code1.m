% Define system parameters
V = [1.06, 1.0, 1.0];
P = [0, 0.65, -3.42];
Q = [0, 0, -0.96];
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
    AP = [];
    AQ = [];

    for i = 1:3
        delta_P = 0;
        delta_Q = 0;

        for j = 1:3
            delta_P = delta_P + V(i) * V(j) * (real(Y(i, j)) * cos(angle(V(i)) - angle(V(j))) + imag(Y(i, j)) * sin(angle(V(i)) - angle(V(j))));
            if i ~= j
                delta_Q = delta_Q + V(i) * V(j) * (real(Y(i, j)) * sin(angle(V(i)) - angle(V(j))) - imag(Y(i, j)) * cos(angle(V(i)) - angle(V(j))));
            end
        end

        if i == 1
            V(i) = V(i) * exp(1j * 0);
        else
            V(i) = sqrt((P(i) - delta_P) * 2 + (Q(i) - delta_Q) * 2) / V(i);
            AP = [AP, P(i) - delta_P];
            AQ = [AQ, Q(i) - delta_Q];
        end
    end

    delta_V = abs(V - V_prev);
    iteration = iteration + 1;
    convergence = [convergence, max(abs(AP)*2 - abs(AQ)*2)];
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

% Create a 3D plot
figure('Position', [100, 100, 800, 600]);
scatter3([1, 0, 0], [0, 0, -1], [0, 0, 0], 200, 'filled', 'MarkerFaceColor', {'red', 'blue', 'green'});
hold on;
plot3([1, 0, 0], [0, 0, -1], [0, 0, 0], '--', 'Color', 'gray');
xlabel('Voltage Magnitude (pu)');
ylabel('Voltage Magnitude (pu)');
zlabel('Voltage Magnitude (pu)');
legend('Bus 1', 'Bus 2', 'Bus 3');
title('Power System Visualization');
grid on;
hold off;