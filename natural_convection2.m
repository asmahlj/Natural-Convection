function h = natural_convection2(A,p, Ts, T_inf)
    % Inputs:
    % D - Diameter or characteristic length (m)
    % Ts - Surface temperature (°C)
    % T_inf - Surrounding air temperature (°C)
    % v - Kinematic viscosity (m^2/s)
    
  [nu, k, Pr] = AirProperties((T_inf+Ts)/2,   [], [],  'nu', 'k', 'Pr');
    
    % Constants
    g = 9.81;  % Gravitational acceleration (m/s^2)
    beta = 1 / ((Ts + T_inf)/2 + 273);  % Coefficient of thermal expansion (1/K)
    % k Thermal conductivity of air (W/m*K)
    
    % Temperature difference
    deltaT = Ts - T_inf;
    
    % Rayleigh number for different configurations
    % Ra = g * beta * deltaT * D^3 / (v^2) * Pr

    
    %% Case 2: Horizontal Plate (facing up or down)
    % Horizontal plate with hot surface facing up
    L_horiz = A/p;  
    
    Ra_horiz_up = g * beta * deltaT * L_horiz^3 / (nu^2) * Pr;
    Nu_horizontal_up = 0.54 * Ra_horiz_up^(1/4);
    h_horizontal_up = Nu_horizontal_up * k / L_horiz;
    
    % Horizontal plate with hot surface facing down
    Nu_horizontal_down = 0.27 * Ra_horiz_up^(1/4);
    h_horizontal_down = Nu_horizontal_down * k / L_horiz;

    

    %% Display Results
    fprintf('Convection Coefficients (h) for Different Shapes:\n');
    fprintf('1. Horizontal Plate (hot surface up): h = %.3f W/m^2K\n', h_horizontal_up);
    fprintf('2. Horizontal Plate (hot surface down): h = %.3f W/m^2K\n', h_horizontal_down);

    
    % Return convection coefficient for all shapes in a struct
    h.horizontal_up = h_horizontal_up;
    h.horizontal_down = h_horizontal_down;

end


