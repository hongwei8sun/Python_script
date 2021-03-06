clear;

parpool('local', str2num(getenv('SLURM_CPUS_PER_TASK')));

%--------------------------------------------------------------
% Define the vertical temperature profile
%--------------------------------------------------------------------

load 'vertical_profile.txt';

Height          = 1000.0*vertical_profile(:,1);   % [m]
P               = vertical_profile(:,2);          % [hPa]
temp            = vertical_profile(:,3);          % [K]
air_num_density = vertical_profile(:,4);          % air number density [molecule/cm3]

concnt1_CO2  = 280.0*10^(-6);
concnt2_CO2  = 560.0*10^(-6);
concnt_H2O   = vertical_profile(:,5)*10^(-6);
concnt_O3    = vertical_profile(:,6)*10^(-6);
concnt_CH4   = vertical_profile(:,7)*10^(-6);
concnt_N2O   = vertical_profile(:,8)*10^(-6);

%---------------------------------------------------------------
% for Q_CO2 
%--------------------------------------------------------------------
file_Q = "Qs_CO2.txt";
fid_Q  = fopen(file_Q);

Q_data = textscan(fid_Q,          ...
    ['%4c' '%22.8f'],               ...
    5000, 'delimiter', '',        ...
    'whitespace', '');
fclose(fid_Q);


Q_look = Q_data{2};

Q_CO2 = Q_look(round(temp));          
% Table: total internal partition sum
% https://hitran.org/docs/iso-meta/

% for Q_CH4 %%%%%%%%%%%%%%%%%%%
file_Q = "Qs_CH4.txt";
fid_Q  = fopen(file_Q);

Q_data = textscan(fid_Q,          ...
    ['%4c' '%22.8f'],               ...
    5000, 'delimiter', '',        ...
    'whitespace', '');
fclose(fid_Q);

Q_look = Q_data{2};

Q_CH4 = Q_look(round(temp)); 

% for Q_H2O %%%%%%%%%%%%%%%%%%%
file_Q = "Qs_H2O.txt";
fid_Q  = fopen(file_Q);

Q_data = textscan(fid_Q,          ...
    ['%4c' '%22.8f'],               ...
    5000, 'delimiter', '',        ...
    'whitespace', '');
fclose(fid_Q);

Q_look = Q_data{2};

Q_H2O = Q_look(round(temp)); 

% for Q_O3 %%%%%%%%%%%%%%%%%%%
file_Q = "Qs_O3.txt";
fid_Q  = fopen(file_Q);

Q_data = textscan(fid_Q,          ...
    ['%4c' '%22.8f'],               ...
    5000, 'delimiter', '',        ...
    'whitespace', '');
fclose(fid_Q);

Q_look = Q_data{2};

Q_O3 = Q_look(round(temp));

% for Q_N2O %%%%%%%%%%%%%%%%%%%
file_Q = "Qs_N2O.txt";
fid_Q  = fopen(file_Q);

Q_data = textscan(fid_Q,          ...
    ['%4c' '%22.8f'],               ...
    5000, 'delimiter', '',        ...
    'whitespace', '');
fclose(fid_Q);

Q_look = Q_data{2};

Q_N2O = Q_look(round(temp));


%-----------------------------------------------------------------
%-----------------------------------------------------------------
g  = 9.8;                 % [m/s2]
R  = 287.058;             % R_specific [J/(kg*K)]
h  = 6.626*10^(-34);      % [J*sec]
c  = 3.0*10^(10);         % [cm/s]
K  = 1.3806*10^(-23);     % [J/K] 
Pi = 3.14;
NA = 6.02214129*10^(23);  % [mol-1]
M_air   = 29;
M_H2O   = 18;
M_CO2   = 44;
%air_rho = air_num_density/NA*M_air  % [g/cm3]

% modify value: $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $
% ***********************************************************
Dv1    = 1;

v1_min = 0.1;
v1_max = 2500.0;

v1     = [v1_min:Dv1:v1_max];                 % corresponding wavelength: 4~100 um

v1_H2O_contimn_B = [v1_min:Dv1:500];
v1_H2O_contimn_T = [v1_min:Dv1:1400];

v1_CO2_contimn_B_a = [v1_min:Dv1:190];
v1_CO2_contimn_T_a = [v1_min:Dv1:450];

v1_CO2_contimn_B_b = [v1_min:Dv1:1150];
v1_CO2_contimn_T_b = [v1_min:Dv1:1800];

X1_CO2   = air_num_density*concnt1_CO2;   % number density [molecule/cm3]
X2_CO2   = air_num_density*concnt2_CO2;    
X_H2O    = air_num_density.*concnt_H2O;    
X_O3     = air_num_density.*concnt_O3;
X_CH4    = air_num_density.*concnt_CH4;
X_N2O    = air_num_density.*concnt_N2O;

rho_H2O   = X_H2O / NA * M_H2O * 10^3;     % [kg/m3]
rho1_CO2  = X1_CO2 / NA * M_CO2 * 10^3;     % [kg/m3]
rho2_CO2  = X2_CO2 / NA * M_CO2 * 10^3;     % [kg/m3]

%------------------------------------------------------------------
% optical depth
%--------------------------------------------------------------------

Sigma1_CO2   = zeros( length(Height) , length(v1) ); % [cm2/molecule]
Sigma2_CO2   = zeros( length(Height) , length(v1) );
Sigma_H2O    = zeros( length(Height) , length(v1) );
Sigma_O3     = zeros( length(Height) , length(v1) );
Sigma_CH4    = zeros( length(Height) , length(v1) );
Sigma_N2O    = zeros( length(Height) , length(v1) );

Kappa_H2O      = zeros( length(Height) , length(v1) ); % [m2/kg]
Kappa_CO2_a    = zeros( length(Height) , length(v1) ); % [m2/kg]
Kappa_CO2_b    = zeros( length(Height) , length(v1) ); % [m2/kg]

delt_height  = Height*0.0;

Dtau1_CO2  = zeros( length(Height) , length(v1) );
Dtau2_CO2  = zeros( length(Height) , length(v1) );
Dtau_H2O   = zeros( length(Height) , length(v1) );
Dtau_O3    = zeros( length(Height) , length(v1) );
Dtau_CH4   = zeros( length(Height) , length(v1) );
Dtau_N2O   = zeros( length(Height) , length(v1) );

Dtau_H2O_contimn   = zeros( length(Height) , length(v1) );
Dtau1_CO2_contimn_a = zeros( length(Height) , length(v1) );
Dtau1_CO2_contimn_b = zeros( length(Height) , length(v1) );
Dtau2_CO2_contimn_a = zeros( length(Height) , length(v1) );
Dtau2_CO2_contimn_b = zeros( length(Height) , length(v1) );

tau1_accumulate_CO2   = zeros( length(Height) , length(v1) );
tau2_accumulate_CO2   = zeros( length(Height) , length(v1) );
tau_accumulate_H2O    = zeros( length(Height) , length(v1) );
tau_accumulate_O3     = zeros( length(Height) , length(v1) );
tau_accumulate_CH4    = zeros( length(Height) , length(v1) );
tau_accumulate_N2O    = zeros( length(Height) , length(v1) );

tau_accumulate_H2O_contimn = zeros( length(Height) , length(v1) );
tau1_accumulate_CO2_contimn_a = zeros( length(Height) , length(v1) );
tau1_accumulate_CO2_contimn_b = zeros( length(Height) , length(v1) );
tau2_accumulate_CO2_contimn_a = zeros( length(Height) , length(v1) );
tau2_accumulate_CO2_contimn_b = zeros( length(Height) , length(v1) );



Dtau1_total   = zeros( length(Height) , length(v1) );
Dtau2_total   = zeros( length(Height) , length(v1) );

tau1_accumulate_total = zeros( length(Height) , length(v1) );
tau2_accumulate_total = zeros( length(Height) , length(v1) );

2
% for H2O and CO2  continuum -----------------------------------------------------
fv_H2O_continm  =  12.17    - 0.051*v1 + 8.32*10^(-5)*v1.^2   - 7.07*10^(-8)*v1.^3   + 2.33*10^(-11)*v1.^4;
fv_CO2_continm_a = -8.853 + 0.028534*v1  - 0.00043194*v1.^2 + 1.4349*10^(-6)*v1.^3  - 1.5539*10^(-9)*v1.^4;
fv_CO2_continm_b = -537.09  + 1.0886*v1   - 0.0007566*v1.^2 + 1.8863*10^(-7)*v1.^3 - 8.2635*10^(-12)*v1.^4;
            
for i=1:length(P)        
    Kappa_H2O(i,:)   = exp(fv_H2O_continm) *(296.0/temp(i))^4.25;  % [m2/kg]
    Kappa_CO2_a(i,:) = exp(fv_CO2_continm_a)*(300.0/temp(i))^1.7;  % [m2/kg]
    Kappa_CO2_b(i,:) = exp(fv_CO2_continm_b)*(300.0/temp(i))^1.7;  % [m2/kg]
end

P_H2O    = concnt_H2O*M_H2O/M_air .* P *100.0; % [Pa]
P1_CO2   = concnt1_CO2*M_CO2/M_air .* P *100.0; % [Pa]
P2_CO2   = concnt2_CO2*M_CO2/M_air .* P *100.0; % [Pa]

for i=2:length(P)

    delt_height(i) = Height(i) - Height(i-1);

    
    H2O_P_rho_B = P_H2O(i-1)/10000.0 * rho_H2O(i-1);
    H2O_P_rho_T = P_H2O(i)/10000.0   * rho_H2O(i);
    CO2_P_rho1_B = P1_CO2(i-1)/10000.0 * rho1_CO2(i-1);
    CO2_P_rho1_T = P1_CO2(i)/10000.0   * rho1_CO2(i);
    CO2_P_rho2_B = P2_CO2(i-1)/10000.0 * rho2_CO2(i-1);
    CO2_P_rho2_T = P2_CO2(i)/10000.0   * rho2_CO2(i);


    H2O_Kappa_P_rho_B  = Kappa_H2O(i-1,:)   * H2O_P_rho_B;
    H2O_Kappa_P_rho_T  = Kappa_H2O(i,:)     * H2O_P_rho_T;

    CO2_Kappa_P_rho1_B_a = Kappa_CO2_a(i-1,:) * CO2_P_rho1_B;
    CO2_Kappa_P_rho1_T_a = Kappa_CO2_a(i,:)   * CO2_P_rho1_T;
    CO2_Kappa_P_rho1_B_b = Kappa_CO2_b(i-1,:) * CO2_P_rho1_B;
    CO2_Kappa_P_rho1_T_b = Kappa_CO2_b(i,:)   * CO2_P_rho1_T;

    CO2_Kappa_P_rho2_B_a = Kappa_CO2_a(i-1,:) * CO2_P_rho2_B;
    CO2_Kappa_P_rho2_T_a = Kappa_CO2_a(i,:)   * CO2_P_rho2_T;
    CO2_Kappa_P_rho2_B_b = Kappa_CO2_b(i-1,:) * CO2_P_rho2_B;
    CO2_Kappa_P_rho2_T_b = Kappa_CO2_b(i,:)   * CO2_P_rho2_T;

    Dtau_H2O_contimn(i,:)  = 0.5*(H2O_Kappa_P_rho_B  + H2O_Kappa_P_rho_T)  * delt_height(i);
    Dtau1_CO2_contimn_a(i,:) = 0.5*(CO2_Kappa_P_rho1_B_a + CO2_Kappa_P_rho1_T_a) * delt_height(i);
    Dtau1_CO2_contimn_b(i,:) = 0.5*(CO2_Kappa_P_rho1_B_b + CO2_Kappa_P_rho1_T_b) * delt_height(i);
    Dtau2_CO2_contimn_a(i,:) = 0.5*(CO2_Kappa_P_rho2_B_a + CO2_Kappa_P_rho2_T_a) * delt_height(i);
    Dtau2_CO2_contimn_b(i,:) = 0.5*(CO2_Kappa_P_rho2_B_b + CO2_Kappa_P_rho2_T_b) * delt_height(i);    


    Dtau_H2O_contimn(i , length(v1_H2O_contimn_T) : length(v1)) = 0.0;
    Dtau_H2O_contimn(i , 1 : length(v1_H2O_contimn_B)         ) = 0.0;

    Dtau1_CO2_contimn_a(i , length(v1_CO2_contimn_T_a) : length(v1)) = 0.0;
    Dtau1_CO2_contimn_a(i , 1 : length(v1_CO2_contimn_B_a)         ) = 0.0;
    Dtau1_CO2_contimn_b(i , length(v1_CO2_contimn_T_b) : length(v1)) = 0.0;
    Dtau1_CO2_contimn_b(i , 1 : length(v1_CO2_contimn_B_b)         ) = 0.0;
    
    Dtau2_CO2_contimn_a(i , length(v1_CO2_contimn_T_a) : length(v1)) = 0.0;
    Dtau2_CO2_contimn_a(i , 1 : length(v1_CO2_contimn_B_a)         ) = 0.0;
    Dtau2_CO2_contimn_b(i , length(v1_CO2_contimn_T_b) : length(v1)) = 0.0;
    Dtau2_CO2_contimn_b(i , 1 : length(v1_CO2_contimn_B_b)         ) = 0.0;

end

3
%--------------------------------------------------------------------------

parfor i=1:length(P)
    Sigma1_CO2(i,:) = absorption_coef_CO2(v1, P(i), temp(i), Q_CO2(i), concnt1_CO2);  % [cm2/molecule]
    Sigma2_CO2(i,:) = absorption_coef_CO2(v1, P(i), temp(i), Q_CO2(i), concnt2_CO2);  % [cm2/molecule]
    Sigma_O3(i,:)   = absorption_coef_O3(v1, P(i), temp(i), Q_O3(i), concnt_O3(i));   % [cm2/molecule]
    Sigma_H2O(i,:)  = absorption_coef_H2O(v1, P(i), temp(i), Q_H2O(i), concnt_H2O(i));
    Sigma_CH4(i,:)  = absorption_coef_CH4(v1, P(i), temp(i), Q_CH4(i), concnt_CH4(i));
    Sigma_N2O(i,:)  = absorption_coef_N2O(v1, P(i), temp(i), Q_N2O(i), concnt_N2O(i));
end

3.5

for i=2:length(P)

    delt_height(i) = Height(i)-Height(i-1);

    Dtau1_CO2(i,:)= 0.5*( Sigma1_CO2(i-1,:)*X1_CO2(i-1) + Sigma1_CO2(i,:) *X1_CO2(i) )* 100.0*delt_height(i);
    Dtau2_CO2(i,:)= 0.5*( Sigma2_CO2(i-1,:)*X2_CO2(i-1) + Sigma2_CO2(i,:) *X2_CO2(i) )* 100.0*delt_height(i);
    Dtau_O3(i,:)  = 0.5*( Sigma_O3(i-1,:)  *X_O3(i-1)   + Sigma_O3(i,:)   *X_O3(i) )  * 100.0*delt_height(i);
    Dtau_H2O(i,:) = 0.5*( Sigma_H2O(i-1,:) *X_H2O(i-1)  + Sigma_H2O(i,:)  *X_H2O(i) ) * 100.0*delt_height(i);
    Dtau_CH4(i,:) = 0.5*( Sigma_CH4(i-1,:) *X_CH4(i-1)  + Sigma_CH4(i,:)  *X_CH4(i) ) * 100.0*delt_height(i);
    Dtau_N2O(i,:) = 0.5*( Sigma_N2O(i-1,:) *X_N2O(i-1)  + Sigma_N2O(i,:)  *X_N2O(i) ) * 100.0*delt_height(i);

    %Dtau1_CO2(i,:)= 0.5*( Sigma1_CO2(i-1,:)+Sigma1_CO2(i,:) )* 0.5*( X1_CO2(i-1)+X1_CO2(i) )* 100.0*delt_height(i);
    %Dtau2_CO2(i,:)= 0.5*( Sigma2_CO2(i-1,:)+Sigma2_CO2(i,:) )* 0.5*( X2_CO2(i-1)+X2_CO2(i) )* 100.0*delt_height(i);
    %Dtau_O3(i,:)  = 0.5*( Sigma_O3(i-1,:)  +Sigma_O3(i,:) )  * 0.5*( X_O3(i-1)  +X_O3(i) )  * 100.0*delt_height(i);
    %Dtau_H2O(i,:) = 0.5*( Sigma_H2O(i-1,:) +Sigma_H2O(i,:) ) * 0.5*( X_H2O(i-1) +X_H2O(i) ) * 100.0*delt_height(i);
    %Dtau_CH4(i,:) = 0.5*( Sigma_CH4(i-1,:) +Sigma_CH4(i,:) ) * 0.5*( X_CH4(i-1) +X_CH4(i) ) * 100.0*delt_height(i);

    Dtau1_total(i,:) =Dtau1_CO2(i,:) +Dtau_O3(i,:) +Dtau_H2O(i,:) +Dtau_CH4(i,:) +Dtau_N2O(i,:) ...
			 +Dtau_H2O_contimn(i,:) +Dtau1_CO2_contimn_a(i,:) +Dtau1_CO2_contimn_b(i,:);
    Dtau2_total(i,:) =Dtau2_CO2(i,:) +Dtau_O3(i,:) +Dtau_H2O(i,:) +Dtau_CH4(i,:) +Dtau_N2O(i,:) ...
			 +Dtau_H2O_contimn(i,:) +Dtau2_CO2_contimn_a(i,:) +Dtau2_CO2_contimn_b(i,:);
    %Dtau1_total(i,:) = Dtau1_CO2(i,:) +Dtau_O3(i,:) +Dtau_H2O(i,:)  +Dtau_CH4(i,:);
    %Dtau2_total(i,:) = Dtau2_CO2(i,:) +Dtau_O3(i,:) +Dtau_H2O(i,:)  +Dtau_CH4(i,:);

    tau1_accumulate_CO2(i,:)         = tau1_accumulate_CO2(i-1,:) + Dtau1_CO2(i,:);
    tau2_accumulate_CO2(i,:)         = tau2_accumulate_CO2(i-1,:) + Dtau2_CO2(i,:);
    tau_accumulate_O3(i,:)           = tau_accumulate_O3(i-1,:)   + Dtau_O3(i,:);
    tau_accumulate_H2O(i,:)          = tau_accumulate_H2O(i-1,:)  + Dtau_H2O(i,:);
    tau_accumulate_CH4(i,:)          = tau_accumulate_CH4(i-1,:)  + Dtau_CH4(i,:);
    tau_accumulate_N2O(i,:)          = tau_accumulate_N2O(i-1,:)  + Dtau_N2O(i,:);
    tau_accumulate_H2O_contimn(i,:)  = tau_accumulate_H2O_contimn(i-1,:)  + Dtau_H2O_contimn(i,:);
    tau1_accumulate_CO2_contimn_a(i,:) = tau1_accumulate_CO2_contimn_a(i-1,:) + Dtau1_CO2_contimn_a(i,:);
    tau1_accumulate_CO2_contimn_b(i,:) = tau1_accumulate_CO2_contimn_b(i-1,:) + Dtau1_CO2_contimn_b(i,:);
    tau2_accumulate_CO2_contimn_a(i,:) = tau2_accumulate_CO2_contimn_a(i-1,:) + Dtau2_CO2_contimn_a(i,:);
    tau2_accumulate_CO2_contimn_b(i,:) = tau2_accumulate_CO2_contimn_b(i-1,:) + Dtau2_CO2_contimn_b(i,:);    

    tau1_accumulate_total(i,:) = tau1_accumulate_CO2(i,:) +tau_accumulate_O3(i,:) +tau_accumulate_H2O(i,:) ...
                        +tau_accumulate_CH4(i,:) +tau_accumulate_N2O(i,:) +tau_accumulate_H2O_contimn(i,:) ...
			+tau1_accumulate_CO2_contimn_a(i,:) +tau1_accumulate_CO2_contimn_b(i,:);
    
    tau2_accumulate_total(i,:) = tau2_accumulate_CO2(i,:) +tau_accumulate_O3(i,:) +tau_accumulate_H2O(i,:) ...
                        +tau_accumulate_CH4(i,:) +tau_accumulate_N2O(i,:) +tau_accumulate_H2O_contimn(i,:) ...
			+tau2_accumulate_CO2_contimn_a(i,:) +tau2_accumulate_CO2_contimn_b(i,:);
                 
end

4
%--------------------------------------------------------------------
% OLR
%--------------------------------------------------------------------

Bv     = Pi*10^4*2.0*h*c*c*v1.*v1.*v1./(exp(h*c*v1./(K*temp(1))) -1.0);

OLR1_CO2 = Bv.*exp( -1.0 * tau1_accumulate_total( length(P),: ) );
OLR2_CO2 = Bv.*exp( -1.0 * tau2_accumulate_total( length(P),: ) );

%--------------------------------------------------------------------

tau1_top  = tau1_accumulate_total(length(P),:);
tau2_top  = tau2_accumulate_total(length(P),:);


Dtau_min = 1.0;
for i=2:length(P)
    i
for j=1:length(v1)
        
    if( Dtau2_total(i,j) > Dtau_min )
            
        Nlev_sub   = floor( Dtau2_total(i,j) / Dtau_min ) + 1;
 
        Dtau1_sub   = Dtau1_total(i,j)/Nlev_sub;
        Dtau2_sub   = Dtau2_total(i,j)/Nlev_sub;

        Dtemp_sub  = ( temp(i)-temp(i-1) )/Nlev_sub;
	
        OLR1_instant = zeros( 1,Nlev_sub );
        OLR2_instant = zeros( 1,Nlev_sub );
            
        parfor i_sub = 1:Nlev_sub

            temp_sub_B = temp(i-1)+Dtemp_sub*(i_sub-1);
            temp_sub_T = temp(i-1)+Dtemp_sub*(i_sub);
            Bv_B      = Pi*10^4*2.0*h*c*c*v1(1,j)*v1(1,j)*v1(1,j)/(exp(h*c*v1(1,j)/(K*temp_sub_B)) -1.0);             
            Bv_T      = Pi*10^4*2.0*h*c*c*v1(1,j)*v1(1,j)*v1(1,j)/(exp(h*c*v1(1,j)/(K*temp_sub_T)) -1.0);             

            tau1_B = tau1_accumulate_total(i-1,j) + Dtau1_sub*(i_sub-1);
            tau1_T = tau1_accumulate_total(i-1,j) + Dtau1_sub*(i_sub);
            tau2_B = tau2_accumulate_total(i-1,j) + Dtau2_sub*(i_sub-1);
            tau2_T = tau2_accumulate_total(i-1,j) + Dtau2_sub*(i_sub);

            fx1_B = Bv_B * exp( tau1_B - tau1_top(j) );
            fx1_T = Bv_T * exp( tau1_T - tau1_top(j) );
            fx2_B = Bv_B * exp( tau2_B - tau2_top(j) );
            fx2_T = Bv_T * exp( tau2_T - tau2_top(j) );

            OLR1_instant(i_sub) = 0.5*( fx1_B + fx1_T ) * Dtau1_sub;
            OLR2_instant(i_sub) = 0.5*( fx2_B + fx2_T ) * Dtau2_sub;

        end

	OLR1_CO2(j) = OLR1_CO2(j) + sum(OLR1_instant);
	OLR2_CO2(j) = OLR2_CO2(j) + sum(OLR2_instant);
        
    else

        temp_B = temp(i-1);
        temp_T = temp(i);
        Bv_B   = Pi*10^4*2.0*h*c*c*v1(1,j)*v1(1,j)*v1(1,j)/(exp(h*c*v1(1,j)/(K*temp_B)) -1.0); 
        Bv_T   = Pi*10^4*2.0*h*c*c*v1(1,j)*v1(1,j)*v1(1,j)/(exp(h*c*v1(1,j)/(K*temp_T)) -1.0); 
    
        tau1_B = tau1_accumulate_total(i-1,j);
        tau1_T = tau1_accumulate_total(i,j);
        tau2_B = tau2_accumulate_total(i-1,j);
        tau2_T = tau2_accumulate_total(i,j);

        fx1_B = Bv_B * exp( tau1_B - tau1_top(j) );
        fx1_T = Bv_T * exp( tau1_T - tau1_top(j) );
        fx2_B = Bv_B * exp( tau2_B - tau2_top(j) );
        fx2_T = Bv_T * exp( tau2_T - tau2_top(j) );

        OLR1_CO2(j) = OLR1_CO2(j) + 0.5*(fx1_B + fx1_T) *Dtau1_total(i,j);
        OLR2_CO2(j) = OLR2_CO2(j) + 0.5*(fx2_B + fx2_T) *Dtau2_total(i,j);
        
    end    
            
end
end

delete(gcp);
% -----------------------------------------------------------------
A = [v1; OLR1_CO2; OLR2_CO2];

fileID = fopen('OLR_CO2.txt','w');
fprintf(fileID,'%15s %15s %15s\n','v1','OLR1_CO2','OLR2_CO2');
fprintf(fileID,'%15.5f %15.5f %15.5f\n',A);
fclose(fileID);

11
sum( (OLR1_CO2-OLR2_CO2) * ( v1(3)-v1(2) ) )

exit;
