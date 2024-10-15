function h = natural_convection(D, Ts, T_inf)
    % Inputs:
    % D - Diameter or characteristic length (m)
    % Ts - Surface temperature (°C)
    % T_inf - Surrounding air temperature (°C)
    % v - Kinematic viscosity (m^2/s)
    
  [nu, k, Pr] = AirProperties((T_inf+Ts)/2,   [], [],  'nu', 'k', 'Pr');
    
    % Constants
    g = 9.81;  % Gravitational acceleration (m/s^2)
    beta = 2 / (Ts + T_inf + 273.15);  % Coefficient of thermal expansion (1/K)
    % k Thermal conductivity of air (W/m*K)
    
    % Temperature difference
    deltaT = Ts - T_inf;
    
    % Rayleigh number for different configurations
    % Ra = g * beta * deltaT * D^3 / (v^2) * Pr

    %% Case 1: Vertical Plate
    L = D;  % Characteristic length for vertical plate is the height (use D as an input)
    Ra_L = g * beta * deltaT * L^3 / (nu^2) * Pr;
    
    Nu_vertical_plate = (0.825+ (0.387*Ra_L^(1/6))/(1+(0.492/Pr)^(9/16))^(8/27))^2;

    h_vertical_plate = Nu_vertical_plate * k / L;



    %% Case 3: Vertical Cylinder
    if D / nu >= 35 * (Pr / nu)^0.25
        % Treat as vertical plate
        Ra_vert_cyl = Ra_L;
        Nu_vert_cyl = Nu_vertical_plate;
    else
        Ra_vert_cyl = g * beta * deltaT * D^3 / (nu^2) * Pr;
        Nu_vert_cyl = (0.825 + 0.387 * Ra_vert_cyl^(1/6) / (1 + (0.492 / Pr)^(9/16))^(8/27))^2;
    end
    h_vert_cyl = Nu_vert_cyl * k / D;

    %% Case 4: Horizontal Cylinder
    Ra_horiz_cyl = g * beta * deltaT * D^3 / (nu^2) * Pr;
    Nu_horiz_cyl = (0.6 + 0.387 * Ra_horiz_cyl^(1/6) / (1 + (0.559 / Pr)^(9/16))^(8/27))^2;
    h_horiz_cyl = Nu_horiz_cyl * k / D;

    %% Case 5: Sphere
    Ra_sphere = g * beta * deltaT * D^3 / (nu^2) * Pr;
    Nu_sphere = 2 + 0.589 * Ra_sphere^(1/4) / (1 + (0.469 / Pr)^(9/16))^(4/9);
    h_sphere = Nu_sphere * k / D;
    
    %% Case 6: Concentric Sphere
    

    %% Display Results
    fprintf('Convection Coefficients (h) for Different Shapes:\n');
    fprintf('1. Vertical Plate: h = %.3f W/m^2K\n', h_vertical_plate);
    fprintf('2. Vertical Cylinder: h = %.3f W/m^2K\n', h_vert_cyl);
    fprintf('3. Horizontal Cylinder: h = %.3f W/m^2K\n', h_horiz_cyl);
    fprintf('4. Sphere: h = %.3f W/m^2K\n', h_sphere);
    
    % Return convection coefficient for all shapes in a struct
    h.vertical_plate = h_vertical_plate;
    h.vertical_cylinder = h_vert_cyl;
    h.horizontal_cylinder = h_horiz_cyl;
    h.sphere = h_sphere;
end


