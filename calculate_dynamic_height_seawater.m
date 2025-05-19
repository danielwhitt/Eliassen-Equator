function dynamic_height = calculate_dynamic_height_seawater(T, S, P)
    % Calculate Dynamic Height from 3D data cubes of Temperature (T), Salinity (S), and Pressure (P)
    % using the "seawater" toolbox with sw_svan for specific volume anomaly.
    %
    % Inputs:
    %   T - 3D array of temperature [°C], dimensions (lon, lat, depth)
    %   S - 3D array of salinity [PSU], dimensions (lon, lat, depth)
    %   P - 1D array of pressure levels [dbar], length equal to the depth dimension
    %
    % Output:
    %   dynamic_height - 3D array of dynamic height [m^2/s^2], relative to the surface (P(1))
    %
    % Dependencies:
    %   Requires the "seawater" toolbox for specific volume anomaly calculation

    % Constants
    g = 9.81; % Acceleration due to gravity [m/s^2]

    % Dimensions
    [nx, ny, nz] = size(T);

    % Preallocate dynamic height array
    dynamic_height = zeros(nx, ny, nz);

    % Loop through each grid point
    for i = 1:nx
        for j = 1:ny
            % Extract vertical profiles of T and S for current grid point
            T_profile = squeeze(T(i, j, :));
            S_profile = squeeze(S(i, j, :));
            P_profile = P(:); % Ensure P is a column vector

            % Compute specific volume anomaly (sw_svan) at this location
            svan = sw_svan(S_profile, T_profile, P_profile);

            % Numerical integration to calculate dynamic height
            % Use cumtrapz to integrate with respect to pressure
            dynamic_height(i, j, :) = cumtrapz(1e4.*P_profile, svan)./g.*100; % dyn-cm
        end
    end
end