clc;
clear;

% Step 1: Define the System
% Load Flow Analysis using Newton-Raphson Method
clc;
clear;

% Step 1: Define the System
% Bus Data: [Bus Type, Voltage Magnitude, Voltage Angle, P, Q]
% Bus Types: 1 = Slack, 2 = PV, 3 = PQ
bus_data = [
    1, 1.06, 0, 0, 0;    % Slack Bus
    2, 1.04, 0, 0.5, 0;  % PV Bus
    3, 1.0, 0, -1.5, -0.8; % PQ Bus
];

% Line Data: [From Bus, To Bus, Resistance (R), Reactance (X)]
line_data = [
    1, 2, 0.02, 0.04;
    1, 3, 0.01, 0.03;
    2, 3, 0.0125, 0.025;
];

% Base Values
baseMVA = 100; % Base power in MVA

% Step 2: Form the Admittance Matrix (Ybus)
n_bus = size(bus_data, 1); % Number of buses
Ybus = zeros(n_bus, n_bus); % Initialize Ybus

for i = 1:size(line_data, 1)
    from = line_data(i, 1);
    to = line_data(i, 2);
    R = line_data(i, 3);
    X = line_data(i, 4);
    Z = R + 1j*X; % Impedance
    Y = 1/Z;      % Admittance
    Ybus(from, to) = Ybus(from, to) - Y;
    Ybus(to, from) = Ybus(to, from) - Y;
    Ybus(from, from) = Ybus(from, from) + Y;
    Ybus(to, to) = Ybus(to, to) + Y;
end

% Step 3: Initialize Variables
V = bus_data(:, 2); % Voltage magnitudes
theta = zeros(n_bus, 1); % Voltage angles (initialized to 0)
P = bus_data(:, 4) / baseMVA; % Real power (pu)
Q = bus_data(:, 5) / baseMVA; % Reactive power (pu)

% Identify PV and PQ buses
pq_buses = find(bus_data(:, 1) == 3); % PQ buses
pv_buses = find(bus_data(:, 1) == 2); % PV buses

% Step 4: Perform Newton-Raphson Iterations
max_iter = 10; % Maximum iterations
tol = 1e-6;    % Tolerance

for iter = 1:max_iter
    % Calculate power mismatches
    [P_calc, Q_calc] = calculate_power(V, theta, Ybus, n_bus);
    dP = P(2:end) - P_calc(2:end); % Exclude slack bus
    dQ = Q(pq_buses) - Q_calc(pq_buses); % Only for PQ buses
    
    % Check for convergence
    if max(abs([dP; dQ])) < tol
        fprintf('Converged in %d iterations.\n', iter);
        break;
    end
    
    % Form the Jacobian matrix
    J = form_jacobian(V, theta, Ybus, n_bus, P_calc, Q_calc, pq_buses, pv_buses);
    
    % Solve for voltage updates
    dX = J \ [dP; dQ];
    dtheta = dX(1:n_bus-1);
    dV = dX(n_bus:end);
    
    % Update voltages
    theta(2:end) = theta(2:end) + dtheta;
    V(pq_buses) = V(pq_buses) + dV; % Only update PQ buses
end

% Step 5: Display Results
disp('Voltage Magnitudes (pu):');
disp(V);
disp('Voltage Angles (radians):');
disp(theta);
disp('Real Power Flows (pu):');
disp(P_calc * baseMVA);
disp('Reactive Power Flows (pu):');
disp(Q_calc * baseMVA);

% Step 6: Calculate Power Flows
function [P_calc, Q_calc] = calculate_power(V, theta, Ybus, n_bus)
    P_calc = zeros(n_bus, 1);
    Q_calc = zeros(n_bus, 1);
    for i = 1:n_bus
        for j = 1:n_bus
            P_calc(i) = P_calc(i) + V(i) * V(j) * (real(Ybus(i,j)) * cos(theta(i)-theta(j)) + imag(Ybus(i,j)) * sin(theta(i)-theta(j)));
            Q_calc(i) = Q_calc(i) + V(i) * V(j) * (real(Ybus(i,j)) * sin(theta(i)-theta(j)) - imag(Ybus(i,j)) * cos(theta(i)-theta(j)));
        end
    end
end

% Step 7: Form the Jacobian Matrix
function J = form_jacobian(V, theta, Ybus, n_bus, P_calc, Q_calc, pq_buses, pv_buses)
    npq = length(pq_buses);
    nvar = (n_bus - 1) + npq; % Total number of variables

    J = zeros(nvar, nvar); % Initialize Jacobian

    % Fill Jacobian matrix
    for i = 2:n_bus
        for j = 2:n_bus
            if i == j
                J(i-1, j-1) = -Q_calc(i) - V(i)^2 * imag(Ybus(i,i)); % dP/dTheta
                if ismember(i, pq_buses)
                    pq_index = find(pq_buses == i); % Get the index in PQ set
                    J(n_bus-1 + pq_index, n_bus-1 + pq_index) = Q_calc(i) / V(i) - V(i) * imag(Ybus(i,i)); % dQ/dV
                end
            else
                J(i-1, j-1) = V(i) * V(j) * (real(Ybus(i,j)) * sin(theta(i)-theta(j)) - imag(Ybus(i,j)) * cos(theta(i)-theta(j))); % dP/dTheta
                if ismember(i, pq_buses)
                    pq_index = find(pq_buses == i);
                    J(n_bus-1 + pq_index, j-1) = -V(i) * V(j) * (real(Ybus(i,j)) * cos(theta(i)-theta(j)) + imag(Ybus(i,j)) * sin(theta(i)-theta(j))); % dQ/dTheta
                end
            end
        end
    end

    % Fill dP/dV and dQ/dV for PQ buses
    for i = 1:npq
        bus = pq_buses(i);
        for j = 2:n_bus
            J(n_bus-1 + i, j-1) = V(bus) * (real(Ybus(bus,j)) * cos(theta(bus)-theta(j)) + imag(Ybus(bus,j)) * sin(theta(bus)-theta(j)));
        end
    end
end
