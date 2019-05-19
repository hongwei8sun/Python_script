function Sigma1 = absorption_coef_CH4(v1, P, T, Q, concnt_CH4) % P:[hPa], T:[K]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c     = 2.99792458*10^(10);  % light speed [cm/s]
c2    = 1.4387770;           % [cm K]
k     = 1.3806488*10^(-16);  % Plank Constant [erg/K]
NA    = 6.02214129*10^(23);  % [mol-1]

% v1    = [2970.0:0.01:2980.0];

P_atm  = P/1013;      % pressure [atm]
P_ref  = 1.0;                % reference pressure 1 [atm]
% T      = 252.4;              % temperature [K]
T_ref  = 296.0;              % [K]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q_ref   = 590.48;         % at 296K
Ma      = 16.0;           % molar mass[g/mol]
P_self = concnt_CH4*Ma/29.0*P_atm;    % Table: partial pressure for a gass [atm]

file_name = ["CH4_313943.par"];

number_of_lines_to_read_per_pass = 110580;  %; 313943;  %
% adjust this number depending on how many lines are in your file


fid = fopen(file_name(1));

HITRAN_data = textscan(fid,                            ...
    ['%2c' '%1f' '%12f' '%10f' '%10f'                  ...
    '%5f' '%5f' '%10f' '%4f' '%8f'                     ...
    '%15c' '%15c' '%15c' '%15c' '%6c'                  ...
    '%12c' '%1c' '%7f' '%7f'],                         ...
    number_of_lines_to_read_per_pass, 'delimiter', '', ...
    'whitespace', '');
fclose(fid);

iso_L      = HITRAN_data{2}; % isotopologue number
nu_L       = HITRAN_data{3}; % line wavenumber [cm^-1]
Sref_L     = HITRAN_data{4}; % line intensity [cm^-1 / (cm^2/molecule)] at Tref
Delt_air   = HITRAN_data{10};% Pressure shift induced by air, referred to p=1 atm [cm-1.atm-1]

gum_air    = HITRAN_data{6}; 
% The air-broadened half width at half maximum (HWHM) [cm?1/atm] at Tref=296K and reference pressure pref=1atm
gum_self   = HITRAN_data{7}; 
%The self-broadened half width at half maximum (HWHM) (cm?1/atm) at Tref=296K and reference pressure pref=1atm
E          = HITRAN_data{8};
%Lower-state energy,[cm-1]
n_air      = HITRAN_data{9}; 
%The coefficient of the temperature dependence of the air-broadened half width

Sigma1  = 0.0.*v1;

for ii=1:number_of_lines_to_read_per_pass

vij = nu_L(ii);      % nu_L
gum = (T_ref/T)^n_air(ii) * ( (gum_air(ii)*(P_atm-P_self)) + gum_self(ii)*P_self );

Sij = Sref_L(ii)*Q_ref/Q *exp(-1.0*c2*E(ii)/T)/exp(-1.0*c2*E(ii)/T_ref) ...
     *(1.0-exp(-1.0*c2*vij/T))/(1.0-exp(-1.0*c2*vij/T_ref)); 

% Sref_L
delt = Delt_air(ii);  % Delt_air

% pressure broaden
fv1   = 1/pi*gum./(gum^2+[v1-(vij+delt*P)].^2);

Sigma1 = Sigma1 + Sij*fv1;
%Sigma1 = Sigma1 + exp(-1.0*c2*E(ii)/T)/exp(-1.0*c2*E(ii)/T_ref);

end  % ii

return