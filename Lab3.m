clear all
clc
D=4e-6; %Calcium Diffusion Coeffecient cm^2/sec
k=1; % Counter for time steps
rstep=0.1*10^-4; %delta r in cm
Ca=linspace(.150,.150,31); % Array representing calcium concentration profile across the cell using 30 slices and beginning concentration of 0.15 micromolar
bc=50; % Boundary concentration of 50 um
j=30;  % Counter for radial positions
dt=0.00001; % Time step

% Main loop for time steps
while k<500 % Run the simulation for 500 time steps
    j=30; % Reset the radial position counter
    Ca(31)=bc; % Set the boundary concentration at the outer edge of the cell, just inside the membrane

     % Loop through radial positions
    while j > 0
        % Calculate change in calcium concentration at each radial position
        if j==1 %center of cell?
            dc(j)=(((4*D)/(rstep^2))*(Ca(2)-Ca(1)))*dt; % Boundary condition at center
            j=j-1; % j becomes 0 to exit the while loop
        else
            % Diffusion equation for calcium concentration for each
            % individual annulus
            dc(j)=((D/rstep^2)*((Ca(j+1)-2*Ca(j)+Ca(j-1))+((Ca(j+1)-Ca(j-1))/(2*(j-1)))))*dt;
            j=(j-1); % move to next annulus getting closer to center
        end
    end
    
    Ca(31)=[]; % Remove the boundary concentration from the array
    Catime(k,:)=Ca+dc; % Store the current calcium concentration, considering changes in concentration
    Ca=Ca+dc; % Update the calcium concentration profile for the next time step
    Ca(31)=bc; % Reupdate the boundary condition
    k=k+1; % proceed to next time step
end

r=linspace(0,3,30); % Radial distance from the center of the cell to 3 um
t=linspace(0,.0049,499); % Time points for which concentration profiles are calculated from 0 to 0.0049 seconds
figure(1)
surf(r,t,Catime) % 3D plot of radium time and concentration
title('3-D Surface Plot of Concentration Profile of Calcium')
xlabel('Distance r (um)')
ylabel('time (seconds)')
zlabel('Concentration (micromolar)')
figure(2)
plot(r,Catime(300,:)) %plot concentration vs radius at t=3 ms
title('Spatial Solution (t=3ms)')
xlabel('Distance r (um)')
ylabel('Concentration (micromolar)')

%%
clear all
clc
D=4e-6; %Calcium Diffusion Coeffecient cm^2/sec
k=1; % Counter for time steps
rstep=0.1*10^-4; %delta r in cm
Ca=linspace(.150,.150,31); % Array representing calcium concentration profile across the cell using 30 slices and beginning concentration of 0.15 micromolar
bc=50; % Boundary concentration of 50 um
j=30;  % Counter for radial positions to start at outer shell
dt=0.00001; % Time step
n = 1;
K_m = 0.500; % Concentration where Vmax/2 occurs (uM)
V_max = 3.2e-7*rstep; % Maximum velocity (uM/s)

% Main loop for time steps
while k<500 % Run the simulation for 500 time steps
    j=30; % Reset the radial position counter
    Ca(31)=bc; % Set the boundary concentration at the outer edge of the cell, just inside the membrane

     % Loop through radial positions
    while j > 0
        % Calculate membrane flux term using Michaelis-Menten kinetics
        pump = V_max * Ca(j) / (K_m + Ca(j));

        % Calculate change in calcium concentration at each radial position
        if j==1 %center of cell
            dc(j)=(((4*D)/(rstep^2))*(Ca(2)-Ca(1)))*dt - pump*dt; % Boundary condition at center
            j=j-1; % j becomes 0 to exit the while loop
        else
            % Diffusion equation for calcium concentration for each
            % individual annulus
            dc(j)=((D/rstep^2)*((Ca(j+1)-2*Ca(j)+Ca(j-1))+((Ca(j+1)-Ca(j-1))/(2*(j-1)))))*dt - pump*dt;
            j=(j-1); % move to next annulus getting closer to center
        end
    end
    
    Ca(31)=[]; % Remove the boundary concentration from the array
    Catime(k,:)=Ca+dc; % Store the current calcium concentration, considering changes in concentration
    Ca=Ca+dc; % Update the calcium concentration profile for the next time step
    Ca(31)=bc; % Reupdate the boundary condition
    k=k+1; % proceed to next time step
end

r=linspace(0,3,30); % Radial distance from the center of the cell to 3 um
t=linspace(0,.0049,499); % Time points for which concentration profiles are calculated from 0 to 0.0049 seconds
figure;
surf(r,t,Catime) % 3D plot of radium time and concentration
title('3-D Surface Plot of Concentration Profile of Calcium')
xlabel('Distance r (um)')
ylabel('time (seconds)')
zlabel('Concentration (micromolar)')
figure;
plot(r,Catime(300,:)) %plot concentration vs radius at t=3 ms
title('Spatial Solution (t=3ms)')
xlabel('Distance r (um)')
ylabel('Concentration (micromolar)')