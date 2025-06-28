clear all; clc; close all;

% MDPE PIPE properties
E =  2e+9;                % Shell material Young's modulus (N/m^2)
rho = 900;                % Density of the shell material (kg/m^3)
etap = 0.06;              % Material Loss factor
a = 0.0845;               % Pipe radius (m) - Shell radius
h = 0.0110;               % Pipe-wall thickness (m)
Ec = E*(1+1j*etap);       % Complex Young's modulus (N/m^2) if the material is lossy

% WATER properties
Bw = 2.25e+9;               % Bulk modulus (N/m^2)
cw = 1500;                  % Speed of sound in water (m/s) - compressional wavespeed
rhow = 1000;                % Density of the water (kg/m^3)

% SOIL A PROPERTIES - SANDY SOIL
Ba = 5.3e+7;                % Bulk modulus (N/m^2)
Ga = 2.00e+7;               % Shear modulus (N/m^2)
rhosa = 2000;               % Density of the soil (kg/m^3)
etaa = 0.0;                 % Bulk and shear loss factor
Bac = Ba*(1+1j*etaa);       % Complex bulk modulus
Gac = Ga*(1+1j*etaa);       % Complex bulk modulus
cd = sqrt((Bac + (4/3)*Gac)/rhosa); % dilatational wavespeed
cs = sqrt(Gac/rhosa);               %Shear wavespeed

% Frequency Vector
dt = 10;
w = 1:dt:600*2*pi;                  % Frequency Vector
kw = w./cw;                         % Wavenumber of water within the pipe (contained in the fluid)

% Other necessary parameters  to estimate wavenumber
Kwater = 2*Bw/a;                    % Dynamic Stiffness of the Water pipe
Kpipe = Ec*h/(a^2) - rho*h*(w.^2);  % Dynamic Stiffness of the pipe-wall
kd2 = (w./cd).^2;                   % Dilatational Wavenumber
kr2 = (w./cs).^2;                   % Shear Wavenumber

% 1) Wavenumber prediction without surrounding soil effects (we use it as
% an initial guess
ka = kw.*sqrt(1 + Kwater./Kpipe);
xa = ka;
Es = 1e-12;
Ea = 100;
xold = xa;
n= 1;                             % Iteration counter

while Ea > Es
    kd1r = sqrt(kd2 - (xa.^2));   % Surrounding medium compressional radial wavenumber
    kr1r = sqrt(kr2 - (xa.^2));   % Surrounding medium shear radial wavenumber
    if abs(real(kd1r)) > abs(imag(kd1r)) & abs(real(kr1r)) > abs(imag(kr1r))
        J0 = besselj(0, kd1r*a);    % Bessel function of the first kind
        Y0 = bessely(0, kd1r*a);    % Bessel function of the second kind
        J1 = besselj(1, kd1r*a);
        Y1 = bessely(1, kd1r*a);
        J0r = besselj(0, kr1r*a);
        Y0r = bessely(0, kr1r*a);
        J1r = besselj(1, kr1r*a);
        Y1r = bessely(1, kr1r*a);

        H = J0 -1j*Y0;              % Hankel function of the second kind, H0 = J0 - iY0
        Hd = -J1 + 1j*Y1;           % Hankel function of the second kind, H0'= -H1 = -J1 + iY1

        Hr = J0r -1j*Y0r;              % Hankel function of the second kind, H0 = J0 - iY0
        Hdr = -J1r + 1j*Y1r;           % Hankel function of the second kind, H0'= -H1 = -J1 + iY1

       K1 = (Bac - (2/3)*Gac).*(kd2./kd1r).*(1 - 2*((xa.^2)./kr2)).*(H./Hd);
       K2 = -2*Gac*kd1r.*(1-(2*((xa.^2)./kr2))).*((H + (Hd./(kd1r*a)))./(-Hd));
       K3 = -4*Gac*kr1r.*((xa.^2)./kr2).*((Hr + (Hdr./(kr1r*a)))./(-Hdr));

    elseif  abs(real(kd1r))<abs(imag(kd1r)) & abs(real(kr1r))>abs(imag(kr1r)) % ii) Shear wave radiates and comp. does not.
        J0 = besselj(0,-kd1r*a);                 % Bessel function of first kind
        Y0 = bessely(0,-kd1r*a);                 % Bessel function of second kind
        J1 = besselj(1,-kd1r*a);                 % Bessel function of first kind
        Y1 = bessely(1,-kd1r*a);                 % Bessel function of second kind
        J0r = besselj(0,kr1r*a);                % Bessel function of first kind
        Y0r = bessely(0,kr1r*a);                % Bessel function of second kind
        J1r = besselj(1,kr1r*a);                % Bessel function of first kind
        Y1r = bessely(1,kr1r*a);                % Bessel function of second kind
        H = J0 - 1j*Y0;                         % Hankel function of the second kind, H0 = J0 - iY0
        Hd = -J1 + 1j*Y1;                       % Hankel function of the second kind, H0' = -H1 = - J1 + iY1
        Hr = J0r - 1j*Y0r;                      % Hankel function of the second kind, H0 = J0 - iY0
        Hdr = -J1r + 1j*Y1r;                    % Hankel function of the second kind, H0' = -H1 = - J1 + iY1
        K1 = (Bac - (2/3)*Gac).*(kd2./-kd1r).*(1 - 2*((xa.^2)./kr2)).*(H./Hd);
        K2 = -2*Gac*-kd1r.*(1-(2*((xa.^2)./kr2))).*((H + (Hd./(-kd1r*a)))./(-Hd));
        K3 = -4*Gac*kr1r.*((xa.^2)./kr2).*((Hr + (Hdr./(kr1r*a)))./(-Hdr));

    elseif abs(real(kd1r))>abs(imag(kd1r)) & abs(real(kr1r))<abs(imag(kr1r)) % iii) Comp. wave radiates and shear wave does not.
        J0 = besselj(0,kd1r*a);                 % Bessel function of first kind
        Y0 = bessely(0,kd1r*a);                 % Bessel function of second kind
        J1 = besselj(1,kd1r*a);                 % Bessel function of first kind
        Y1 = bessely(1,kd1r*a);                 % Bessel function of second kind
        J0r = besselj(0,-kr1r*a);                % Bessel function of first kind
        Y0r = bessely(0,-kr1r*a);                % Bessel function of second kind
        J1r = besselj(1,-kr1r*a);                % Bessel function of first kind
        Y1r = bessely(1,-kr1r*a);                % Bessel function of second kind
        H = J0 - 1j*Y0;                         % Hankel function of the second kind, H0 = J0 - iY0
        Hd = -J1 + 1j*Y1;                       % Hankel function of the second kind, H0' = -H1 = - J1 + iY1
        Hr = J0r - 1j*Y0r;                      % Hankel function of the second kind, H0 = J0 - iY0
        Hdr = -J1r + 1j*Y1r;                    % Hankel function of the second kind, H0' = -H1 = - J1 + iY1
        K1 = (Bac - (2/3)*Gac).*(kd2./kd1r).*(1 - 2*((xa.^2)./kr2)).*(H./Hd);
        K2 = -2*Gac*kd1r.*(1-(2*((xa.^2)./kr2))).*((H + (Hd./(kd1r*a)))./(-Hd));
        K3 = -4*Gac*-kr1r.*((xa.^2)./kr2).*((Hr + (Hdr./(-kr1r*a)))./(-Hdr));

    else                                                                      % iv) Both shear and comp. waves do not radiate.
        J0 = besselj(0,-kd1r*a);                % Bessel function of first kind
        Y0 = bessely(0,-kd1r*a);                % Bessel function of second kind
        J1 = besselj(1,-kd1r*a);                % Bessel function of first kind
        Y1 = bessely(1,-kd1r*a);                % Bessel function of second kind
        J0r = besselj(0,-kr1r*a);               % Bessel function of first kind
        Y0r = bessely(0,-kr1r*a);               % Bessel function of second kind
        J1r = besselj(1,-kr1r*a);               % Bessel function of first kind
        Y1r = bessely(1,-kr1r*a);               % Bessel function of second kind
        H = J0 - 1j*Y0;                         % Hankel function of the second kind, H0 = J0 - iY0
        Hd = -J1 + 1j*Y1;                       % Hankel function of the second kind, H0' = -H1 = - J1 + iY1
        Hr = J0r - 1j*Y0r;                      % Hankel function of the second kind, H0 = J0 - iY0
        Hdr = -J1r + 1j*Y1r;                    % Hankel function of the second kind, H0' = -H1 = - J1 + iY1
        K1 = (Bac - (2/3)*Gac).*(kd2./-kd1r).*(1 - 2*((xa.^2)./kr2)).*(H./Hd);
        K2 = -2*Gac*-kd1r.*(1-(2*((xa.^2)./kr2))).*((H + (Hd./(-kd1r*a)))./(-Hd));
        K3 = -4*Gac*-kr1r.*((xa.^2)./kr2).*((Hr + (Hdr./(-kr1r*a)))./(-Hdr));
    end
%==================== Function associated to the implicit equation of
%wavenumber prediction ==================%
xa = kw.*sqrt(1 + (Kwater./(Kpipe+K1+K2+K3)));
Ea = (norm((xa - xold)./xa))* 100;
xold = xa;
n = n+1;
end
%=============== Wave attenuation ===============%
wattexa = -20*imag(xa)./log(10);

%=============== Wavespeed analysis =============%
cxa = (w./real(xa));
% Dispersion curves for the s=1 wave in a MDPE pipe buried in SANDY soil - Wavenumber, Attenuation and normalized wavespeed
%=============== Plot the Result ================%
figure(1)
subplot(1,2,1)
plot(w/2/pi,real(xa));hold on
ylim([0 18]);
xlabel('Frequency [Hz]');
ylabel('Real part of wavenumber (1/m)');
subplot(1,2,2)
plot(w/2/pi,wattexa);hold on
ylim([0 20]);
xlabel('Frequency [Hz]');
ylabel('Attenuation (dB/m)');

%% ======================= Monte Carlo Simulation ========================= %%
% Set the random seed
rng(1000);
% Define the mean and standard deviation of the parameters
mu_a = a; sigma_a = 0.01 * mu_a;
mu_h = h; sigma_h = 0.01 * mu_h;
mu_Ec = Ec; sigma_Ec = 0.01 * mu_Ec;
mu_rho = rho; sigma_rho = 0.01 * mu_rho;
mu_Bac = Bac; sigma_Bac = 0.01 * mu_Bac;
mu_Gac = Gac; sigma_Gac = 0.01 * mu_Gac;
mu_rhosa = rhosa; sigma_rhosa = 0.01 * mu_rhosa;

tic
% Number of repetitions or simulations
rep = 1000;
% Generate Gaussian Distributed Samples for each parameter
norm_a = normrnd(mu_a, sigma_a, [rep, 1]);
norm_h = normrnd(mu_h, sigma_h, [rep, 1]);
norm_Ec = normrnd(mu_Ec, sigma_Ec, [rep, 1]);
norm_rho = normrnd(mu_rho, sigma_rho, [rep, 1]);
norm_Bac = normrnd(mu_Bac, sigma_Bac, [rep, 1]);
norm_Gac = normrnd(mu_Gac, sigma_Gac, [rep, 1]);
norm_rhosa = normrnd(mu_rhosa, sigma_rhosa, [rep, 1]);

% Initialize cell array to store output (wave speeds and attenuation) for each realization
rc = cell(rep, 1);
rwa = cell(rep, 1);

% Loop through the simulations
for i = 1:rep
    cd = sqrt((norm_Bac(i)+(4/3)*norm_Gac(i))/norm_rhosa(i)); % Dilatational wavespeed
    cs = sqrt(norm_Gac(i)/norm_rhosa(i));
    Kwater =  2*Bw/norm_a(i);                   % Dynamic stiffness of the water pipe (pressure/disp)
    Kpipe = norm_Ec(i)*norm_h(i)/(norm_a(i).^2) - norm_rho(i).*norm_h(i).*(w.^2);  % Dynamic stiffneess of the pipe-wall
    kd2 = (w./cd).^2;                            % Compressional (dilatational) wavenumber
    kr2 = (w./cs).^2;
    ka = kw.*sqrt(1 + Kwater./Kpipe); % in air
    xa  = ka;                                    % Initial guess based on ka from Korteweg Eq
    Es = 1e-12;                                  % Tolerance
    Ea = 100;                                     % Randomly large relative approximate error
    xold = xa;
    n = 1;                                       % Iteration counter

    while Ea > Es
        kd1r = sqrt(kd2 - (xa.^2));                   % Surrounding medium compressional radial wavenumbers
        kr1r = sqrt(kr2 - (xa.^2));                   % Surrounding medium shear radial wavenumbers
        if abs(real(kd1r))>abs(imag(kd1r)) & abs(real(kr1r))>abs(imag(kr1r)) %  i) Both shear and compressional waves radiate.
            J0 = besselj(0,kd1r*norm_a(i));                 % Bessel function of first kind
            Y0 = bessely(0,kd1r*norm_a(i));                 % Bessel function of second kind
            J1 = besselj(1,kd1r*norm_a(i));                 % Bessel function of first kind
            Y1 = bessely(1,kd1r*norm_a(i));                 % Bessel function of second kind
            J0r = besselj(0,kr1r*norm_a(i));                % Bessel function of first kind
            Y0r = bessely(0,kr1r*norm_a(i));                % Bessel function of second kind
            J1r = besselj(1,kr1r*norm_a(i));                % Bessel function of first kind
            Y1r = bessely(1,kr1r*norm_a(i));                % Bessel function of second kind
            H = J0 - 1j*Y0;                         % Hankel function of the second kind, H0 = J0 - iY0
            Hd = -J1 + 1j*Y1;                       % Hankel function of the second kind, H0' = -H1 = - J1 + iY1
            Hr = J0r - 1j*Y0r;                      % Hankel function of the second kind, H0 = J0 - iY0
            Hdr = -J1r + 1j*Y1r;                    % Hankel function of the second kind, H0' = -H1 = - J1 + iY1
            K1 = (norm_Bac(i) - (2/3)*norm_Gac(i)).*(kd2./kd1r).*(1 - 2*((xa.^2)./kr2)).*(H./Hd);
            K2 = -2*norm_Gac(i)*kd1r.*(1-(2*((xa.^2)./kr2))).*((H + (Hd./(kd1r*norm_a(i))))./(-Hd));
            K3 = -4*norm_Gac(i)*kr1r.*((xa.^2)./kr2).*((Hr + (Hdr./(kr1r*norm_a(i))))./(-Hdr));

        elseif abs(real(kd1r))<abs(imag(kd1r)) & abs(real(kr1r))>abs(imag(kr1r)) % ii) Shear wave radiates and comp. does not.
            J0 = besselj(0,-kd1r*norm_a(i));                 % Bessel function of first kind
            Y0 = bessely(0,-kd1r*norm_a(i));                 % Bessel function of second kind
            J1 = besselj(1,-kd1r*norm_a(i));                 % Bessel function of first kind
            Y1 = bessely(1,-kd1r*norm_a(i));                 % Bessel function of second kind
            J0r = besselj(0,kr1r*norm_a(i));                % Bessel function of first kind
            Y0r = bessely(0,kr1r*norm_a(i));                % Bessel function of second kind
            J1r = besselj(1,kr1r*norm_a(i));                % Bessel function of first kind
            Y1r = bessely(1,kr1r*norm_a(i));                % Bessel function of second kind
            H = J0 - 1j*Y0;                         % Hankel function of the second kind, H0 = J0 - iY0
            Hd = -J1 + 1j*Y1;                       % Hankel function of the second kind, H0' = -H1 = - J1 + iY1
            Hr = J0r - 1j*Y0r;                      % Hankel function of the second kind, H0 = J0 - iY0
            Hdr = -J1r + 1j*Y1r;                    % Hankel function of the second kind, H0' = -H1 = - J1 + iY1
            K1 = (norm_Bac(i) - (2/3)*norm_Gac(i)).*(kd2./-kd1r).*(1 - 2*((xa.^2)./kr2)).*(H./Hd);
            K2 = -2*norm_Gac(i)*-kd1r.*(1-(2*((xa.^2)./kr2))).*((H + (Hd./(-kd1r*norm_a(i))))./(-Hd));
            K3 = -4*norm_Gac(i)*kr1r.*((xa.^2)./kr2).*((Hr + (Hdr./(kr1r*norm_a(i))))./(-Hdr));

        elseif abs(real(kd1r))>abs(imag(kd1r)) & abs(real(kr1r))<abs(imag(kr1r)) % iii) Comp. wave radiates and shear wave does not.
            J0 = besselj(0,kd1r*norm_a(i));                 % Bessel function of first kind
            Y0 = bessely(0,kd1r*norm_a(i));                 % Bessel function of second kind
            J1 = besselj(1,kd1r*norm_a(i));                 % Bessel function of first kind
            Y1 = bessely(1,kd1r*norm_a(i));                 % Bessel function of second kind
            J0r = besselj(0,-kr1r*norm_a(i));                % Bessel function of first kind
            Y0r = bessely(0,-kr1r*norm_a(i));                % Bessel function of second kind
            J1r = besselj(1,-kr1r*norm_a(i));                % Bessel function of first kind
            Y1r = bessely(1,-kr1r*norm_a(i));                % Bessel function of second kind
            H = J0 - 1j*Y0;                         % Hankel function of the second kind, H0 = J0 - iY0
            Hd = -J1 + 1j*Y1;                       % Hankel function of the second kind, H0' = -H1 = - J1 + iY1
            Hr = J0r - 1j*Y0r;                      % Hankel function of the second kind, H0 = J0 - iY0
            Hdr = -J1r + 1j*Y1r;                    % Hankel function of the second kind, H0' = -H1 = - J1 + iY1
            K1 = (norm_Bac(i) - (2/3)*norm_Gac(i)).*(kd2./kd1r).*(1 - 2*((xa.^2)./kr2)).*(H./Hd);
            K2 = -2*norm_Gac(i)*kd1r.*(1-(2*((xa.^2)./kr2))).*((H + (Hd./(kd1r*norm_a(i))))./(-Hd));
            K3 = -4*norm_Gac(i)*-kr1r.*((xa.^2)./kr2).*((Hr + (Hdr./(-kr1r*norm_a(i))))./(-Hdr));

        else                                                                      % iv) Both shear and comp. waves do not radiate.
            J0 = besselj(0,-kd1r*norm_a(i));                % Bessel function of first kind
            Y0 = bessely(0,-kd1r*norm_a(i));                % Bessel function of second kind
            J1 = besselj(1,-kd1r*norm_a(i));                % Bessel function of first kind
            Y1 = bessely(1,-kd1r*norm_a(i));                % Bessel function of second kind
            J0r = besselj(0,-kr1r*norm_a(i));               % Bessel function of first kind
            Y0r = bessely(0,-kr1r*norm_a(i));               % Bessel function of second kind
            J1r = besselj(1,-kr1r*norm_a(i));               % Bessel function of first kind
            Y1r = bessely(1,-kr1r*norm_a(i));               % Bessel function of second kind
            H = J0 - 1j*Y0;                         % Hankel function of the second kind, H0 = J0 - iY0
            Hd = -J1 + 1j*Y1;                       % Hankel function of the second kind, H0' = -H1 = - J1 + iY1
            Hr = J0r - 1j*Y0r;                      % Hankel function of the second kind, H0 = J0 - iY0
            Hdr = -J1r + 1j*Y1r;                    % Hankel function of the second kind, H0' = -H1 = - J1 + iY1
            K1 = (norm_Bac(i) - (2/3)*norm_Gac(i)).*(kd2./-kd1r).*(1 - 2*((xa.^2)./kr2)).*(H./Hd);
            K2 = -2*norm_Gac(i)*-kd1r.*(1-(2*((xa.^2)./kr2))).*((H + (Hd./(-kd1r*norm_a(i))))./(-Hd));
            K3 = -4*norm_Gac(i)*-kr1r.*((xa.^2)./kr2).*((Hr + (Hdr./(-kr1r*norm_a(i))))./(-Hdr));
        end

    % 1) Function associated to the implicit equation of wavenumber prediction
    xa = kw.*sqrt(1 + (Kwater./(Kpipe+K1+K2+K3)));

    Ea = (norm((xa - xold)./xa))* 100;
    xold = xa;
    n = n+1;
  end

    % 2) Wave attenuation
    auxwa = -20*imag(xa)./log(10);

    % 3) Wavespeed analysis
    auxc =  (w./real(xa));

    rc{i} = auxc;
    rwa{i} = auxwa;
end

rcm = cell2mat(rc);
rwam = cell2mat(rwa);

%% Plot the simulation
left = 100; bottom = 100; width = 500; height = 400;

figure
b = fill([w/2/pi fliplr(w/2/pi)],[max(rcm) fliplr(min(rcm))],'g');
set(b,'EdgeColor', [.7 .7 .8], 'FaceColor', [.7 .7 .8]);hold on
plot(w/2/pi, cxa,'b--', 'LineWidth', 2.0);
ylim([300 450]);
xlabel('Frequency [Hz]');
ylabel('Wavespeed (m/s)');
set(gca, 'FontName', 'Times', 'Fontsize', 14);
ax = gca;
ax.XColor = 'k';
ax.YColor = 'k';
myfiguresize = [left, bottom, width, height];

figure
set(gcf, 'Position', myfiguresize); 
b = fill([w/2/pi fliplr(w/2/pi)],[max(rwam) fliplr(min(rwam))],'g');
set(b,'EdgeColor', [.7 .7 .8], 'FaceColor', [.7 .7 .8]);hold on
plot(w/2/pi,wattexa,'b--', 'LineWidth', 2.0);
ylim([0 30]);
xlabel('Frequency [Hz]');
ylabel('Attenuation (dB/m)');
set(gca, 'FontName', 'Times', 'Fontsize', 14);
ax = gca;
ax.XColor = 'k';
ax.YColor = 'k';
myfiguresize = [left, bottom, width, height]; 

toc
% ============================== Mean  ============================ %
mean_wavespeed = mean(rcm, 1);
mean_attenuation = mean(rwam, 1);

% ====================== Standard Deviation ======================= %
std_wavespeed = std(rcm, 0, 1);
std_attenuation = std(rwam, 0, 1);

% ===================== Percentile ================================ %

% Wavespeed
percentile_5_wavespeed = prctile(rcm, 5, 1);
percentile_95_wavespeed = prctile(rcm, 95, 1);
% Attenuation
percentile_5_attenuation = prctile(rwam, 5, 1);
percentile_95_attenuation = prctile(rwam, 95, 1);

freq_index = 50; 

figure;
subplot(1,2,1); 
histogram(rcm(:, freq_index)); % Histogram of wavespeed at freq_index 
title(['Wave Speed Distribution at ', num2str(w(freq_index)/(2*pi)), ' Hz']);
xlabel('Wave Speed (m/s)');
ylabel('Frequency Count');

subplot(1,2,2);
histogram(rwam(:, freq_index)); % Histogram of attenuation at freq_index 
title(['Attenuation Distribution at ', num2str(w(freq_index)/(2*pi)), ' Hz']);
xlabel('Attenuation (dB/m)');
ylabel('Frequency Count');

% =========================== Confidence Interval ===================== %
N = rep; % Number of realizations
Z_score = 1.96; % Z-score for a 95% confidence interval

% For Wave Speed
mean_ws_at_freq = mean_wavespeed(freq_index);
std_ws_at_freq = std_wavespeed(freq_index);
sem_ws = std_ws_at_freq / sqrt(N); % Standard Error of the Mean
ci_ws_lower = mean_ws_at_freq - Z_score * sem_ws;
ci_ws_upper = mean_ws_at_freq + Z_score * sem_ws;

% For Attenuation
mean_att_at_freq = mean_attenuation(freq_index);
std_att_at_freq = std_attenuation(freq_index);
sem_att = std_att_at_freq / sqrt(N);
ci_att_lower = mean_att_at_freq - Z_score * sem_att;
ci_att_upper = mean_att_at_freq + Z_score * sem_att;

fprintf('At %.2f Hz:\n', w(freq_index)/(2*pi));
fprintf('  95%% Confidence Interval for Mean Wave Speed: [%.3f, %.3f] m/s\n', ci_ws_lower, ci_ws_upper);
fprintf('  95%% Confidence Interval for Mean Wave Attenuation: [%.3f, %.3f] dB/m\n', ci_att_lower, ci_att_upper); 

figure('Name', 'Wave Speed with 90% Confidence Interval'); 
set(gcf, 'Position', myfiguresize); 

% Plot 5th and 95th percentile envelope for wave speed
fill([w/2/pi fliplr(w/2/pi)],[percentile_95_wavespeed fliplr(percentile_5_wavespeed)],'c','FaceAlpha',0.3);
hold on;
plot(w/2/pi, mean_wavespeed, 'r-', 'LineWidth', 2); % Plot the mean
xlabel('Frequency [Hz]');
ylabel('Wave Speed (m/s)');
title('Wave Speed with 90% Confidence Interval');
legend('90% Confidence Interval', 'Mean Wave Speed', 'Location', 'best'); 
grid on;

figure('Name', 'Attenuation with 90% Confidence Interval'); 
set(gcf, 'Position', myfiguresize); 
% Plot 5th and 95th percentile envelope for wave attenuation
fill([w/2/pi fliplr(w/2/pi)], [percentile_95_attenuation fliplr(percentile_5_attenuation)], 'm', 'FaceAlpha', 0.3);
hold on;
plot(w/2/pi, mean_attenuation, 'r-', 'LineWidth', 2);
xlabel('Frequency [Hz]');
ylabel('Attenuation [dB/m]'); 
title('Wave Attenuation with 90% Confidence Interval');
legend('90% Confidence Interval', 'Mean Wave Attenuation', 'Location', 'best'); 
grid on;

% ============= Sensitivity Analysis (Pearson Correlation) ============= %
output_wavespeed_freq = rcm(:, freq_index);
% Input parameters for correlation
inputs_for_correlation = [
    norm_a, ...
    norm_h, ...
    real(norm_Ec), imag(norm_Ec), ...
    norm_rho, ...
    real(norm_Bac), imag(norm_Bac), ...
    real(norm_Gac), imag(norm_Gac), ...
    norm_rhosa ...
];
% Calculate the Pearson Correlation Matrix
corr_matrix = corrcoef([inputs_for_correlation, output_wavespeed_freq]);

% Extract correlations with the output (last column of the last row, excluding the last element itself)
correlation_coefficients = corr_matrix(end, 1:end-1);

% Define input parameter names for clear display
input_param_names = {
    'Pipe Radius (a)', 'Pipe Thickness (h)', 'Ec (Real)', 'Ec (Imag)', ...
    'Pipe Density (rho)', 'Bac (Real)', 'Bac (Imag)', 'Gac (Real)', 'Gac (Imag)', 'Soil Density (rhosa)'
};
fprintf('\nSensitivity Analysis (Pearson Correlation to Wave Speed) at %.2f Hz:\n', w(freq_index)/(2*pi));
% Display sorted results
[sorted_corr, sorted_indices] = sort(abs(correlation_coefficients), 'descend');
for k = 1:length(sorted_indices)
    idx = sorted_indices(k);
    fprintf('  %s: %.4f\n', input_param_names{idx}, correlation_coefficients(idx));
end

%% ======================================= Probability Density Function ====================================== %%
% For Wave Speed PDF at freq_index
mu_ws_freq = mean_wavespeed(freq_index);
sigma_ws_freq = std_wavespeed(freq_index);
x_ws = linspace(min(rcm(:, freq_index)), max(rcm(:, freq_index)), 100); 
pdf_ws = normpdf(x_ws, mu_ws_freq, sigma_ws_freq);

figure;
plot(x_ws, pdf_ws, 'LineWidth', 2)
xlabel('Wavespeed (m/s)')
ylabel('Probability Density')
title(['Normal Distribution PDF for Wave Speed at ', num2str(w(freq_index)/(2*pi)), ' Hz'])
grid on;

% For Wave Attenuation PDF at freq_index
mu_att_freq = mean_attenuation(freq_index);
sigma_att_freq = std_attenuation(freq_index);
x_att = linspace(min(rwam(:, freq_index)), max(rwam(:, freq_index)), 100); 
pdf_att = normpdf(x_att, mu_att_freq, sigma_att_freq);

figure;
plot(x_att, pdf_att, 'LineWidth', 2)
xlabel('Wave Attenuation (dB/m)')
ylabel('Probability Density')
title(['Normal Distribution PDF for Wave Attenuation at ', num2str(w(freq_index)/(2*pi)), ' Hz'])
grid on;


%% ======================= Monte Carlo Convergence Analysis ========================= %%

% Maximum number of repetitions 
max_rep = rep; 

% Initialize arrays to store the evolving mean and standard deviation
mean_wavespeed_convergence = zeros(1, max_rep);
std_wavespeed_convergence = zeros(1, max_rep);
mean_attenuation_convergence = zeros(1, max_rep);
std_attenuation_convergence = zeros(1, max_rep);


% Loop through increasing number of simulations to observe convergence
for num_sims = 1:max_rep
    % Mean and standard deviation for the current number of simulations for the selected frequency index
    mean_wavespeed_convergence(num_sims) = mean(rcm(1:num_sims, freq_index));
    std_wavespeed_convergence(num_sims) = std(rcm(1:num_sims, freq_index));
    mean_attenuation_convergence(num_sims) = mean(rwam(1:num_sims, freq_index));
    std_attenuation_convergence(num_sims) = std(rwam(1:num_sims, freq_index));
end

% Plotting Convergence
figure('Name', 'Monte Carlo Convergence');
left = 100; bottom = 100; width = 1000; height = 800;
set(gcf, 'Position', [left, bottom, width, height]);

subplot(2,2,1);
plot(1:max_rep, mean_wavespeed_convergence, 'b-', 'LineWidth', 1.5);
xlabel('Number of Simulations');
ylabel('Mean Wave Speed (m/s)');
title(['Convergence of Mean Wave Speed at ', num2str(w(freq_index)/(2*pi)), ' Hz']);
grid on;

subplot(2,2,2);
plot(1:max_rep, std_wavespeed_convergence, 'r-', 'LineWidth', 1.5);
xlabel('Number of Simulations');
ylabel('Standard Deviation of Wave Speed (m/s)');
title(['Convergence of Std Dev of Wave Speed at ', num2str(w(freq_index)/(2*pi)), ' Hz']);
grid on;

subplot(2,2,3);
plot(1:max_rep, mean_attenuation_convergence, 'b-', 'LineWidth', 1.5);
xlabel('Number of Simulations');
ylabel('Mean Attenuation (dB/m)');
title(['Convergence of Mean Attenuation at ', num2str(w(freq_index)/(2*pi)), ' Hz']);
grid on;

subplot(2,2,4);
plot(1:max_rep, std_attenuation_convergence, 'r-', 'LineWidth', 1.5);
xlabel('Number of Simulations');
ylabel('Standard Deviation of Attenuation (dB/m)');
title(['Convergence of Std Dev of Attenuation at ', num2str(w(freq_index)/(2*pi)), ' Hz']);
grid on;

% Error analysis
final_mean_ws = mean_wavespeed_convergence(end);
final_std_ws = std_wavespeed_convergence(end);
final_mean_att = mean_attenuation_convergence(end);
final_std_att = std_attenuation_convergence(end);

error_mean_ws = abs(mean_wavespeed_convergence - final_mean_ws);
error_std_ws = abs(std_wavespeed_convergence - final_std_ws);
error_mean_att = abs(mean_attenuation_convergence - final_mean_att);
error_std_att = abs(std_attenuation_convergence - final_std_att);

figure('Name', 'Monte Carlo Error');
set(gcf, 'Position', [left, bottom, width, height]);

subplot(2,2,1);
semilogy(1:max_rep, error_mean_ws, 'b-', 'LineWidth', 1.5);
xlabel('Number of Simulations');
ylabel('Absolute Error in Mean Wave Speed');
title(['Error in Mean Wave Speed Convergence at ', num2str(w(freq_index)/(2*pi)), ' Hz']);
grid on;

subplot(2,2,2);
semilogy(1:max_rep, error_std_ws, 'r-', 'LineWidth', 1.5);
xlabel('Number of Simulations');
ylabel('Absolute Error in Std Dev of Wave Speed');
title(['Error in Std Dev of Wave Speed Convergence at ', num2str(w(freq_index)/(2*pi)), ' Hz']);
grid on;

subplot(2,2,3);
semilogy(1:max_rep, error_mean_att, 'b-', 'LineWidth', 1.5);
xlabel('Number of Simulations');
ylabel('Absolute Error in Mean Attenuation');
title(['Error in Mean Attenuation Convergence at ', num2str(w(freq_index)/(2*pi)), ' Hz']);
grid on;

subplot(2,2,4);
semilogy(1:max_rep, error_std_att, 'r-', 'LineWidth', 1.5);
xlabel('Number of Simulations');
ylabel('Absolute Error in Std Dev of Attenuation');
title(['Error in Std Dev of Attenuation Convergence at ', num2str(w(freq_index)/(2*pi)), ' Hz']);
grid on;
