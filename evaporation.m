% Set the year and calculate the total number of days in the year
year = 2023;
t = yeardays(year);
HOUR = t * 24;

% Load regridded data for initial temperature and MERRA-2 data
load('Dec_AvgTemp_Regridded.mat', 'avgTempGrid');
load(['MERRA2_regridded_data_', num2str(year), '.mat']);

% Determine the field name for initial temperature based on the year
initial_year_field = sprintf('Year%d', year - 1);
initial_temp_data = avgTempGrid.(initial_year_field);

% Load the area shapefile for the MENA region
area = ['MENA_area_', num2str(year), '.shp'];

% Read and process the shapefile data
Table1 = readgeotable(area);
T1 = geotable2table(Table1, ["Latitude", "Longitude"]);
T1 = sortrows(T1, "ID");
Data1 = T1.total_sum;         % Read the area attribute (in m^2 per cell)
Data1(Data1 == 0) = nan;      % Set zero areas to NaN
Ax = reshape(Data1, [216, 103]);
Area1 = rot90(Ax);            % Rotate the area matrix

% Initialize cell array for daily evaporation values and matrices for monthly and yearly values
dailyEvaporation = cell(103, 216);
monthlyEvaporation = nan(103, 216, 12);
yearlyEvaporation = nan(103, 216);

% Loop through each grid cell
for i=1:103
    for j=1:216
        if Area1(i, j) > 0

            % Extract relevant data from the regridded MERRA-2 dataset
            U = MERRA2_SPEEDLMLq_1(i, j, :);         % Wind speed (m/s)
            Tair = MERRA2_TLMLq_1(i, j, :);          % Air temperature (K)
            Rs = MERRA2_SWGDNq_1(i, j, :);           % Incoming radiation (W/m^2)
            Qa = MERRA2_QLMLq_1(i, j, :);            % Specific humidity (kg/kg)
            PS = MERRA2_PSq_1(i, j, :);              % Air pressure (Pa)

            % Reshape data for interpolation
            U_1 = U(:);
            Tair_1 = Tair(:);
            Rs_1 = Rs(:);
            Qa_1 = Qa(:);
            PS_1 = PS(:);

            % Interpolate data to a finer time resolution (half-hourly)
            U_data = interp1(1:HOUR, U_1, 1:0.5:HOUR);
            Tair_data = interp1(1:HOUR, Tair_1, 1:0.5:HOUR);
            Rs_data = interp1(1:HOUR, Rs_1, 1:0.5:HOUR);
            Qa_data = interp1(1:HOUR, Qa_1, 1:0.5:HOUR);
            PS_data = interp1(1:HOUR, PS_1, 1:0.5:HOUR);

            % Calculate vapor pressure (Pa)
            e_v = 28.97 / 18 * Qa_data .* PS_data;

            %%% Calculation of evaporation  

            % Constants for the evaporation model

            Z = 0:0.04:3;               % Depth intervals (m)
            ro = 1000;                  % Water density (kg/m^3)
            c = 4200;                   % Heat capacity of water (J/kgK)
            alpha = 0.143 * 10^-6;      % Molecular thermal diffusivity of water (m^2/s)
            e_w = 0.95;                 % Water surface emissivity (-)
            albedo = 0.05;              % Water surface albedo (-)
            Beta = 0.3;                 % Absorption coefficient at the water surface (-)
            eta = 0.02;                 % Radiation attenuation coef. (-)
            Da = 2.5 * 10^-5;           % Vapor diffusion coefficient in air (m^2/s)
            L = 2450000;                % Latent heat of water vaporization (J/kg)
            k_air = 0.0257;             % Thermal conductivity of the air (W/mK)

            % Linearization of saturated vapor concentration: Cs=a1*Tw+b1 (278-298 K)
            a1 = 0.0008;
            b1 = -0.2169;

            e_a = 0.8;                  % Emissivity of the air (-)
            Ste = 5.67 * 10^-8;         % Stefan-Boltzmann coefficient (W/m^2K^4)
            T_bot = 299;                % Bottom temperature (K)
            q_b = 0;                    % Heat flux from soil layer beneath (W/m^2)
            dH = 12;
            T_soil = 15 + 273;

            Size_Z = size(Z);
            I = Size_Z(2);              % Number of nodes
            day = t -1;                 % Number of days
            DT = 1;                     % Time step difference
            N = 86400 / DT;             % Number of time steps per day
            dt = 1 * DT;                % Time step (s)

            mm=1;
            for dd=1:day                % The loop for the entire period

                % Initial condition of the water body
                if dd==1
                    T(1:I,1) = initial_temp_data(i, j);  % Set initial temp.
                else
                    T(1:I,1) = T_day(1:I,dd-1);    % Use previous day's temperature profile
                end
               
                count=0;
                for n=1:N-1
                    t1=(n)*dt;
                    count=count+1;
                    if count==1800/dt      % From half hour to second
                        mm=mm+1;
                        count=1;
                    end

                    % To interpolate climate data from half-hourly to every time step
                    U=(U_data(mm+1)-U_data(mm))/(1800/dt)*count+U_data(mm);
                    delta_m=0.0126/U;                % Thickness of air boundary layer (m)
                    ha=k_air/delta_m;                % Sensible heat flux coefficient (W/m^2K)
                    w=1.2*10^-3*U;                   % Friction velocity of the air (m/s)
                    k_fc=0.4;                        % von Karman's constant
                    k_star=6*U^(-1.84)*1;
                    P0=1;                                     % Prandtl number

                    Rs=(Rs_data(mm+1)-Rs_data(mm))/(1800/dt)*count+Rs_data(mm);
                    PS=(PS_data(mm+1)-PS_data(mm))/(1800/dt)*count+PS_data(mm);
                    Qa=(Qa_data(mm+1)-Qa_data(mm))/(1800/dt)*count+Qa_data(mm);
                    Ta=((Tair_data(mm+1)-Tair_data(mm))/(1800/dt)*count+Tair_data(mm));
                    if Ta<273
                        Ta=273;
                    end

                    % Calculate air vapor pressure and vapor concentration
                    e0 = 101.325 * exp(13.3185 * (1 - 373.15 / Ta) - 1.9760 * (1 - 373.15 / Ta)^2 - 0.6445 * (1 - 373.15 / Ta)^3 - 0.1299 * (1 - 373.15 / Ta)^4);
                    Ca = (28.97 / 18 * Qa * PS) / (462 * Ta);

                    % Solve temperature profile except surface and bottom
                    for x=2:I-1
                        z=Z(x);
                        dz1=Z(x)-Z(x-1);     % Distance between neighbouring nodes
                        dz2=Z(x+1)-Z(x);

                        ro(x)=(1-1.9549*10^-5*(abs(T(x,n)-277))^1.68)*10^3;    % Water density 
                        N2=9.8/ro(x)*abs((ro(x)-ro(x-1))/dz1);  
                        if N2==0
                            Ri=-1/20;
                        else
                            Ri=(-1+(1+40*N2*(k_fc*z/(w*exp(-k_star*z)))^2)^0.5)/20;
                        end
                        D_tdif=k_fc*w*z/P0*exp(-k_star*z)*(1+37*Ri^2)^-1*(1); 
                        if D_tdif>(5*10^-4)
                            D_tdif=5*10^-4;
                        end
                        Diffius(x)=alpha+D_tdif;
                        T(x,n+1)=T(x,n)+dt*((alpha+D_tdif)/(dz1/2+dz2/2)*((T(x+1,n)-T(x,n))/dz2+(T(x-1,n)-T(x,n))/dz1)+((1+37*Ri^2)^-1)*(k_fc*w/P0*exp(-k_star*z)-k_fc*w*z*k_star/P0*exp(-k_star*z))*(T(x+1,n)-T(x,n))/dz2+((1-Beta)*(1-albedo)*Rs/(ro(x)*c)*(0.032*0.237*exp(eta*exp(-eta*z)))));
                    end

                    ro(1)=(1-1.9549*10^-5*(abs(T(1,n)-277))^1.68)*10^3;
                    ro(I)=(1-1.9549*10^-5*(abs(T(I,n)-277))^1.68)*10^3;

                    % Solve for the surface
                    dz=Z(2)-Z(1);
                    T(1,n+1)=1/(Da*L/delta_m*a1+ha+4*Ste*e_w*Ta^3+(ro(1)*c*alpha)/dz)*((ro(1)*c*alpha)/dz*T(2,n+1)+Beta*(1-albedo)*Rs+Ste*(3*e_w+e_a)*Ta^4+ha*Ta+Da*L/delta_m*(Ca-b1));
                    % Solve for the bottom
                    zo=Z(I);
                    dz=Z(I)-Z(I-1);
                    Rs_b=(1-Beta)*(1-albedo)*Rs*(exp(-eta*zo));
                    T(I,n+1)=1/((ro(I)*c*alpha)/dz+dH)*((ro(I)*c*alpha)/dz*T(I-1,n+1)+dH*T_soil+Rs_b);

                    Time1(n)=(t1-dt)/(3600);     % Save time
                    T_surf(n)=T(1,n);            % Save water surface temperature

                    % Calculate evaporation: LE (W/m^2)
                    % Vapor pressure at water surface
                    e0=101.325*exp(13.3185*(1-373.15./T_surf(n))-1.9760*(1-373.15./T_surf(n)).^2-.6445*(1-373.15./T_surf(n)).^3-.1299*(1-373.15./T_surf(n)).^4);   
                    % Moisture concentration at water surface
                    Cs=e0*1000/(462*T_surf(n));                                                                                                                    
                    LE(n)=Da*L/delta_m*(Cs-Ca)/28.4;
                    % Calculate sensible heat flux: H (W/m^2)
                    H(n)=ha*(T_surf(n)-Ta);

                    if n==(N-1)
                        T_surf(n+1)=T(1,n+1);
                        Time1(n+1)=t1/3600;
                    end

                end

                % Implement the mixing in the reservoir
                for x=1:I
                    T_day(x,dd)=mean(T(x,1:N));
                    T_day1(x,dd)=mean(T(x,1:N));
                end

                [Tmax,idx] = max(T_day(1:I,dd));

                if (idx>2)&&(idx<I)&&(T_day(idx,dd)>T_day(1,dd)+0.05)&&(T_day(idx,dd)>(T_day(I,dd)+0.1))

                    rr=0;
                    for cc=-6:6
                        rr=rr+1;
                        IDx=round(idx/2)+cc;

                        if (IDx>1)&&(IDx<idx)

                            Temp=T_day(1:I,dd);
                            Plc=find(Temp(IDx)<Temp);
                            S=size(Plc);

                            Q1 = trapz(Z(1:IDx),(T_day(IDx,dd)-T_day(1:IDx,dd)));
                            Q2 = trapz(Z(IDx:Plc(S(1))),(T_day(IDx:Plc(S(1)),dd)-T_day(IDx,dd)));
                            Dife(rr)=abs(Q2-Q1);
                            Plc_Dife(rr)=IDx;
                            PLC(rr)= Plc(S(1));
                        else
                            Dife(rr)=20;
                            Plc_Dife(rr)=IDx;
                            PLC(rr)= 0;
                        end

                    end

                    [MIN_D,idx_min] = min(Dife);
                    T_day(1:(PLC(idx_min)),dd)=T_day(Plc_Dife(idx_min),dd);

                elseif (idx>2)&&(T_day(idx,dd)>T_day(1,dd)+0.05)&&((T_day(idx,dd)-0.1)<=T_day(I,dd))
                    Q3 = trapz(Z(1:I),(T_day(1:I,dd)-T_day(1,dd)));
                    base=T_day(1,dd)+Q3/Z(I);
                    for x=1:I
                        T_day(x,dd)=base;
                    end

                end
                Time_day(dd)=dd;
                LE_day1(dd)=mean(LE(1000:(N-1000)));
                dd
            end
            % Store daily evaporation values
            dailyEvaporation{i, j} = LE_day1;
        else
            % Store NaN for cells with no area
            dailyEvaporation{i, j} = NaN;     
        end
        %%% End of evaporation calculations
    end
end

% Save the daily, monthly, and yearly evaporation values
save(['xdailyEvaporation_MENA_4m_', num2str(year), '.mat'], 'dailyEvaporation');
% Save('monthlyEvaporation_MENA_2016.mat', 'monthlyEvaporation');
% Save('yearlyEvaporation_MENA_2016.mat', 'yearlyEvaporation');

