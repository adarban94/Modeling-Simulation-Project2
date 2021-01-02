clc;format long
%% Defining the require constants of catalyts and reactor dimensions

L = 2.5;                         % [m]                              Length of reactor.
dt = 0.0256;                     % [m]                              Diameter of reactor.
dp = 0.0082;                     % [m]                              Diameter of catalyts particle.
epsilon = 0.48;                  % [m^3/m^3]                        Porosity.
Density_bed = 50;                % [kg/m^3]                         Density of reactor bed.
Phis = 1;                        % [unitless]                       Sphericity factor for spheres catalyts particles.
as = 6*(1-epsilon)/(Phis*dp);    % [m^2/m^3]                        External surface to particle volume ratio.

%% Defining the require constants of operation conditions

R = 8.2057e-5;                         % [(atm.m3)/(mol*K)]                      Gas constant.
Pt = 1;                                % [atm]                            Pressure of reactor bed.
Tb = 450+273.15;                       % [K]                              Temperature of reactor coolant.
T0 = 300+273.15;                       % [K]                              Temperature of inlet reactor.
Rep = 1400;                            % [unitless]                       Reynolds number.
Flowin = 4;                            % [Nm^3/h]                         Inlet volume flowrate.
u0=Flowin/((pi/4)*dt^2);               % [m/h]                            Superficial velocity.
y_Air_in = 0.98;                       % [%mol]                           Mole frac of inlet Air.
y_N2_in = 0.79*y_Air_in;               % [%mol]                           Mole frac of inlet Nitrogen.               ### Ref: Che-galicia 2018
y_C2H6_in = 0.02;                      % [%mol]                           Mole frac of inlet Ethane.(1-40 %)         ### Ref: Che-galicia 2018
y_C2H4_in = 0;                         % [%mol]                           Mole frac of inlet Ethene.
y_O2_in = 0.21*y_Air_in;               % [%mol]                           Mole frac of inlet Oxygen.                 ### Ref: Che-galicia 2018
y_CO2_in = 0;                          % [%mol]                           Mole frac of inlet Carbon dioxid.
y_CO_in = 0;                           % [%mol]                           Mole frac of inlet Carbon monoxid.
y_H2O_in = 0;                          % [%mol]                           Mole frac of inlet Water.
y = [y_C2H6_in          ...
    y_C2H4_in y_O2_in  ...
    y_CO2_in  y_CO_in  ...
    y_H2O_in  y_N2_in     ];          % [%mol]                           Mole frac list of total componets [C2H6 C2H4 O2 CO2 CO H2O N2]
C0_C2H6  =  ((Pt*y_C2H6_in)/(R*T0));   % [mol/m^3]                        Inlet concentration of C2H6
C0_C2H4  =  ((Pt*y_C2H4_in)/(R*T0));   % [mol/m^3]                        Inlet concentration of C2H4
C0_O2    =  ((Pt*y_O2_in)/(R*T0))  ;   % [mol/m^3]                        Inlet concentration of O2
C0_CO2   =  ((Pt*y_CO2_in)/(R*T0)) ;   % [mol/m^3]                        Inlet concentration of CO2
C0_CO    =  ((Pt*y_CO_in)/(R*T0))  ;   % [mol/m^3]                        Inlet concentration of CO
C0_H2O   =  ((Pt*y_H2O_in)/(R*T0)) ;   % [mol/m^3]                        Inlet concentration of H2O
C0_N2    =  ((Pt*y_N2_in)/(R*T0))  ;   % [mol/m^3]                        Inlet concentration of N2
C0 = [C0_C2H6          ...
    C0_C2H4   C0_O2  ...
    C0_CO2    C0_CO  ...
    C0_H2O    C0_N2     ];           % [mol/m^3]                        Inlet concentration

%% Defining the require constants

Deffr = 32;                      % [m^2/h]                          Effective mass transfer
%                                  coefficient in radius direction.
Deffz = 53;                      % [m^2/h]                          Effective mass transfer coefficient
%                                  in horizontal axis direction ### Ref: Che-galicia 2018.

kg = 576;                        % [m^3/(m^2*h)]                    Surface mass transfer coefficient.
hg = 928.8;                      % [kJ/(m^2*h*K)]                   Surface heat transfer coefficient.

keffz = 10;  %Assumption         %[kJ/(m*h*K)]                      Effective thermal conductivity.
%                                  in the radius direction.
keffr = 9.72;                    %[kJ/(m*h*K)]                      Effective thermal conductivity.
%                                  in the radius direction, ESTIMATED BY BOUNDARY LAYER APPROX ### Ref: Che-galicia 2018.

hw = 1051.2;                     % [kJ/(m^2*h*K)]                   Wall heat transfer coefficient, ESTIMATED BY BOUNDARY LAYER APPROX ### Ref: Che-galicia 2018.

%% Defining the require constants of components properties
%
% Units
%
% Mw: [g/mol]           Tc: [K]       Pc: [atm]
% cp_R: [unitless depend on R]        deltaS0: [J/(mol*K)]      deltaH0:[kJ/mol]
R=8.314e-3; % (kJ/mol*K)

C2H6 = struct('Mw',30.07,   'Tc',305.406,   'Pc',4880109/101325,...
    'cp_R',[1.131,0.019225,-0.000005561,0,1500] ,'deltaS0',-5.27e01 ,'deltaH0',-4.80e01);
C2H4 = struct('Mw',28.054,  'Tc',282.3438,  'Pc',5045427/101325,...
    'cp_R',[1.424,0.014394,-0.000004392,0,1500] ,'deltaS0',-4.34e01 ,'deltaH0',-1.48e02);
O2   = struct('Mw',31.998,  'Tc',154.645,   'Pc',5043213/101325,...
    'cp_R',[3.639,0.000506,0,-22700,2000]       ,'deltaS0',-5.59e01 ,'deltaH0',-6.02e01);
CO2  = struct('Mw',44.009,  'Tc',304.1548,  'Pc',7380862/101325,...
    'cp_R',[5.457,0.001045,0,-115700,2000]      ,'deltaS0',-5.66e01 ,'deltaH0',-8.38e01);
CO   = struct('Mw',28.01,   'Tc',134.18,    'Pc',3710046/101325,...
    'cp_R',[3.376,0.000557,0,-3100,2500]        ,'deltaS0',-8.66e01 ,'deltaH0',-4.09e01);
H2O  = struct('Mw',18.015,  'Tc',647.1081,  'Pc',22072227/101325,...
    'cp_R',[3.47,0.00145,0,12100,2000]          ,'deltaS0',-5.27e01 ,'deltaH0',-8.63e01);
N2   = struct('Mw',28.014,  'Tc',126.2069,  'Pc',3398154/101325,...
    'cp_R',[3.28,0.000593,0,4000,2000]          ,'deltaS0',0        ,'deltaH0',0);

Components = [C2H6 C2H4 O2 CO2 CO H2O N2];       % List of components

%% Defining the require constants for kinetic of reactions
%
% Units
%
% Aprime(A'): [mmol/(g*h)]   EnergyA: [kJ/mol]      m: [unitless]
% component coefficients  vcoffrxn: [unitless]   deltaHstd: [kJ/mol]
% component order list: [C2H6 C2H4 O2 CO2 CO H2O N2]

RxnKinetic = struct('Aprime',[4.95 1.35 1.76 2.61 2.16],...
    'EnergyA', [7.55e01 5.24e02 1.43e02 1.1e02 8.8e01],...
    'm', [1 5.45e-02 1.07 1.71e-01 5.38e-01],...
    'deltaHstd', [-111.43 -1443.15 -860 -1331.81 -760],...
    'vcoffrxn', [-1 1 -0.5 0 0 1 0; -1 0 -3.5 2 0 3 0;...
    -1 0 -2.5 0 2 3 0; 0 -1 -3 2 0 2 0; 0 -1 -2 0 2 2 0]);

%% Calculation

%===Interior points and coefficients matrix -------------------------------

Nz = 3; % No. of interior point in z direction.
Nr = 1 ; % No. of interior point in r direction.
zmin = 0; zmax = L;
rmin = 0; rmax = dt/2;
z_nodes = [0,sort(Roots_of_Jacobi_Polynomial(0,0,Nz))',1] ;  % Roots of Jacobi polynomial with (a,b==0) in z direction.
% z_nodes = (zmax-zmin)*z_nodes+zmin;
r_nodes = [0,sort(Roots_of_Jacobi_Polynomial(0,0,Nr))',1] ;  % Roots of Jacobi polynomial with (a,b==0) in r direction.
% r_nodes = (rmax-rmin)*r_nodes+rmin;
syms z
Lz = sym(ones(numel(z_nodes),1));
for i=1:numel(z_nodes)
    for j=1:numel(z_nodes)
        if j~=i
            Lz(i,1) = (z-z_nodes(j))/(z_nodes(i)-z_nodes(j))*Lz(i,1);
            % Lz is Lagrange interpolation polynomial in z direction.
        end
    end
end
syms r
Lr = sym(ones(numel(r_nodes),1));
for i=1:numel(r_nodes)
    for j=1:numel(r_nodes)
        if j~=i
            Lr(i,1) = (r-r_nodes(j))/(r_nodes(i)-r_nodes(j))*Lr(i,1);
            % Lr is Lagrange interpolation polynomial in r direction.
        end
    end
end
Lz_prime = diff(Lz);      % First drivative of Lagrange polynomial in z direction.
Lz_Zegond = diff(Lz,2);   % Second drivative of Lagrange polynomial in z direction.
Lr_prime = diff(Lr);      % First drivative of Lagrange polynomial in r direction.
Lr_Zegond = diff(Lr,2);   % Second drivative of Lagrange polynomial in r direction.
Az = zeros(numel(z_nodes));
Bz = zeros(numel(z_nodes));
for i = 1:numel(z_nodes)
    for j = 1:numel(z_nodes)
        Az(i,j) = double(subs(Lz_prime(j),z_nodes(i)));
        Bz(i,j) = double(subs(Lz_Zegond(j),z_nodes(i)));
    end
end
Ar = zeros(numel(r_nodes));
Br = zeros(numel(r_nodes));
for i = 1:numel(r_nodes)
    for j = 1:numel(r_nodes)
        Ar(i,j) = double(subs(Lr_prime(j),r_nodes(i)));
        Br(i,j) = double(subs(Lr_Zegond(j),r_nodes(i)));
    end
end

% load Orthogonal_Matrix.mat
%===Initial guess ---------------------------------------------------------

Nz = length(z_nodes)-2;  % Declare the number of Interior nodes for BC
Nr = length(r_nodes)-2;  % Declare the number of Interior nodes for BC
R_nodes = r_nodes(2:end-1);

Initial_Guess_C_C2H6        =  ones(Nz,Nr)*C0_C2H6;
Initial_Guess_C_C2H4        =  ones(Nz,Nr)*C0_C2H4;
Initial_Guess_C_O2          =  ones(Nz,Nr)*C0_O2;
Initial_Guess_C_CO2         =  ones(Nz,Nr)*C0_CO2;
Initial_Guess_C_CO          =  ones(Nz,Nr)*C0_CO;
Initial_Guess_C_H2O         =  ones(Nz,Nr)*C0_H2O;
Initial_Guess_C_N2          =  ones(Nz,Nr)*C0_N2;
Initial_Guess_Rof           =  ones(Nz,Nr)*0.6159;                   % [kg/(m^3)] Aspen Hysys at inlet condition
Initial_Guess_Cs_C2H6       =  Initial_Guess_C_C2H6;
Initial_Guess_Cs_C2H4       =  Initial_Guess_C_C2H4;
Initial_Guess_Cs_O2         =  Initial_Guess_C_O2  ;
Initial_Guess_Cs_CO2        =  Initial_Guess_C_CO2 ;
Initial_Guess_Cs_CO         =  Initial_Guess_C_CO  ;
Initial_Guess_Cs_H2O        =  Initial_Guess_C_H2O ;
Initial_Guess_Cpf           =  ones(Nz,Nr)*1.084 ;                   % [kJ/(kg.C)] Aspen Hysys at inlet condition
Initial_Guess_T             =  ones(Nz,Nr)*T0;
Initial_Guess_Ts            =  ones(Nz,Nr)*T0;

Initial_Guess=[reshape(Initial_Guess_C_C2H6,1,Nz*Nr)  ,  reshape(Initial_Guess_C_C2H4,1,Nz*Nr)          ,...
    reshape(Initial_Guess_C_O2,1,Nz*Nr)    ,  reshape(Initial_Guess_C_CO2,1,Nz*Nr)           ,...
    reshape(Initial_Guess_C_CO,1,Nz*Nr)    ,  reshape(Initial_Guess_C_H2O,1,Nz*Nr)           ,...
    reshape(Initial_Guess_C_N2,1,Nz*Nr)    ,  reshape(Initial_Guess_Rof,1,Nz*Nr)             ,...
    reshape(Initial_Guess_Cs_C2H6,1,Nz*Nr) ,  reshape(Initial_Guess_Cs_C2H4,1,Nz*Nr)         ,...
    reshape(Initial_Guess_Cs_O2,1,Nz*Nr)   ,  reshape(Initial_Guess_Cs_CO2,1,Nz*Nr)          ,...
    reshape(Initial_Guess_Cs_CO,1,Nz*Nr)   ,  reshape(Initial_Guess_Cs_H2O,1,Nz*Nr)          ,...
    reshape(Initial_Guess_Cpf,1,Nz*Nr)     ,  reshape(Initial_Guess_T,1,Nz*Nr)               ,...
    reshape(Initial_Guess_Ts,1,Nz*Nr)     ];

%===Solver ----------------------------------------------------------------

Option = optimoptions('fsolve','Algorithm','levenberg-marquardt',...
    'Display','iter','MaxFunctionEvaluations',20000000,'MaxIterations',2000000,...
    'FunctionTolerance',1e-20,'StepTolerance',1e-10);

X=fsolve(@(x) Equations(x,Az,Bz,Ar,Br,Nz,Nr,u0,C0,Pt,T0,hw,Tb,epsilon,Density_bed,...
    Flowin,Deffz,Deffr,keffz,keffr,R_nodes,kg,hg,as,Components,RxnKinetic,R),Initial_Guess,Option);

C_C2H6_hat   = sym(zeros(1)) ;       C_C2H4_hat  = sym(zeros(1))  ;
C_O2_hat     = sym(zeros(1)) ;       C_CO2_hat   = sym(zeros(1))  ;
C_CO_hat     = sym(zeros(1)) ;       C_H2O_hat   = sym(zeros(1))  ;
C_N2_hat     = sym(zeros(1)) ;       Rof_hat     = sym(zeros(1))  ;
Cs_C2H6_hat  = sym(zeros(1)) ;       Cs_C2H4_hat = sym(zeros(1))  ;
Cs_O2_hat    = sym(zeros(1)) ;       Cs_CO2_hat  = sym(zeros(1))  ;
Cs_CO_hat    = sym(zeros(1)) ;       Cs_H2O_hat  = sym(zeros(1))  ;
Cpf_hat      = sym(zeros(1)) ;       T_hat       = sym(zeros(1))  ;
Ts_hat       = sym(zeros(1)) ;

C_C_C2H6  = zeros(numel(z_nodes),numel(r_nodes))    ;   C_C_C2H4  = zeros(numel(z_nodes),numel(r_nodes))  ;
C_C_O2    = zeros(numel(z_nodes),numel(r_nodes))    ;   C_C_CO2   = zeros(numel(z_nodes),numel(r_nodes))  ;
C_C_CO    = zeros(numel(z_nodes),numel(r_nodes))    ;   C_C_H2O   = zeros(numel(z_nodes),numel(r_nodes))  ;
C_C_N2    = zeros(numel(z_nodes),numel(r_nodes))    ;   C_Rof     = zeros(numel(z_nodes),numel(r_nodes))  ;
C_Cs_C2H6 = zeros(numel(z_nodes),numel(r_nodes))    ;   C_Cs_C2H4 = zeros(numel(z_nodes),numel(r_nodes))  ;
C_Cs_O2   = zeros(numel(z_nodes),numel(r_nodes))    ;   C_Cs_CO2  = zeros(numel(z_nodes),numel(r_nodes))  ;
C_Cs_CO   = zeros(numel(z_nodes),numel(r_nodes))    ;   C_Cs_H2O  = zeros(numel(z_nodes),numel(r_nodes))  ;
C_Cpf     = zeros(numel(z_nodes),numel(r_nodes))    ;   C_T       = zeros(numel(z_nodes),numel(r_nodes))  ;
C_Ts      = zeros(numel(z_nodes),numel(r_nodes))    ;

C_C_C2H6(2:end-1,2:end-1)  = reshape(X(1:Nz*Nr),Nz,Nr)              ;   C_C_C2H4(2:end-1,2:end-1)  = reshape(X(Nz*Nr+1:2*Nz*Nr),Nz,Nr)      ;
C_C_O2(2:end-1,2:end-1)    = reshape(X(2*Nz*Nr+1:3*Nz*Nr),Nz,Nr)    ;   C_C_CO2(2:end-1,2:end-1)   = reshape(X(3*Nz*Nr+1:4*Nz*Nr),Nz,Nr)    ;
C_C_CO(2:end-1,2:end-1)    = reshape(X(4*Nz*Nr+1:5*Nz*Nr),Nz,Nr)    ;   C_C_H2O(2:end-1,2:end-1)   = reshape(X(5*Nz*Nr+1:6*Nz*Nr),Nz,Nr)    ;
C_C_N2(2:end-1,2:end-1)    = reshape(X(6*Nz*Nr+1:7*Nz*Nr),Nz,Nr)    ;   C_Rof(2:end-1,2:end-1)     = reshape(X(7*Nz*Nr+1:8*Nz*Nr),Nz,Nr)    ;
C_Cs_C2H6(2:end-1,2:end-1) = reshape(X(8*Nz*Nr+1:9*Nz*Nr),Nz,Nr)    ;   C_Cs_C2H4(2:end-1,2:end-1) = reshape(X(9*Nz*Nr+1:10*Nz*Nr),Nz,Nr)   ;
C_Cs_O2(2:end-1,2:end-1)   = reshape(X(10*Nz*Nr+1:11*Nz*Nr),Nz,Nr)  ;   C_Cs_CO2(2:end-1,2:end-1)  = reshape(X(11*Nz*Nr+1:12*Nz*Nr),Nz,Nr)  ;
C_Cs_CO(2:end-1,2:end-1)   = reshape(X(12*Nz*Nr+1:13*Nz*Nr),Nz,Nr)  ;   C_Cs_H2O(2:end-1,2:end-1)  = reshape(X(13*Nz*Nr+1:14*Nz*Nr),Nz,Nr)  ;
C_Cpf(2:end-1,2:end-1)     = reshape(X(14*Nz*Nr+1:15*Nz*Nr),Nz,Nr)  ;   C_T(2:end-1,2:end-1)       = reshape(X(15*Nz*Nr+1:16*Nz*Nr),Nz,Nr)  ;
C_Ts(2:end-1,2:end-1)      = reshape(X(16*Nz*Nr+1:17*Nz*Nr),Nz,Nr)  ;

%Boundry condition

%Initial guess for solving nonlieanr boundary alg equations
C_solid_In = C0(1:6);                 % mole concentration of components in order: [C2H6 C2H4 O2 CO2 CO H2O]
C_solid_Out= C0(1:6);
C_Cpf_In   = 1.084;                   % [kJ/(kg.C)] Aspen Hysys at inlet condition
C_Cpf_Out  = 1.084;                   % [kJ/(kg.C)] Aspen Hysys at inlet condition
C_Rof_In   = 0.6159;                  % [kg/(m^3)] Aspen Hysys at inlet condition
C_Rof_Out  = 0.6159;                  % [kg/(m^3)] Aspen Hysys at inlet condition
Ts0        = T0;

for k = 2:numel(r_nodes)-1
    %z=0
    C_C_C2H6(1,k)  = (((u0*C0(1))+(epsilon*Deffz*Az(1,2:end-1)*C_C_C2H6(2:end-1,k))-(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,2:end-1)*C_C_C2H6(2:end-1,k)))/((u0)+(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,1))-(epsilon*Deffz*Az(1,1))));
    C_C_C2H4(1,k)  = (((u0*C0(2))+(epsilon*Deffz*Az(1,2:end-1)*C_C_C2H4(2:end-1,k))-(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,2:end-1)*C_C_C2H4(2:end-1,k)))/((u0)+(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,1))-(epsilon*Deffz*Az(1,1))));
    C_C_O2(1,k)    = (((u0*C0(3))+(epsilon*Deffz*Az(1,2:end-1)*C_C_O2(2:end-1,k))-(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,2:end-1)*C_C_O2(2:end-1,k)))/((u0)+(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,1))-(epsilon*Deffz*Az(1,1))));
    C_C_CO2(1,k)   = (((u0*C0(4))+(epsilon*Deffz*Az(1,2:end-1)*C_C_CO2(2:end-1,k))-(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,2:end-1)*C_C_CO2(2:end-1,k)))/((u0)+(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,1))-(epsilon*Deffz*Az(1,1))));
    C_C_CO(1,k)    = (((u0*C0(5))+(epsilon*Deffz*Az(1,2:end-1)*C_C_CO(2:end-1,k))-(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,2:end-1)*C_C_CO(2:end-1,k)))/((u0)+(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,1))-(epsilon*Deffz*Az(1,1))));
    C_C_H2O(1,k)   = (((u0*C0(6))+(epsilon*Deffz*Az(1,2:end-1)*C_C_H2O(2:end-1,k))-(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,2:end-1)*C_C_H2O(2:end-1,k)))/((u0)+(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,1))-(epsilon*Deffz*Az(1,1))));
    C_C_N2(1,k)    = (((u0*C0(7))+(epsilon*Deffz*Az(1,2:end-1)*C_C_N2(2:end-1,k))-(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,2:end-1)*C_C_N2(2:end-1,k)))/((u0)+(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,1))-(epsilon*Deffz*Az(1,1))));
    C_gas_In       = [C_C_C2H6(1,k) C_C_C2H4(1,k) C_C_O2(1,k) C_C_CO2(1,k) C_C_CO(1,k) C_C_H2O(1,k) C_C_N2(1,k)];         % mole fraction of components in order: [C2H6 C2H4 O2 CO2 CO H2O N2]
    C_gas_Out      = [C_C_C2H6(end,k) C_C_C2H4(end,k) C_C_O2(end,k) C_C_CO2(end,k) C_C_CO(end,k) C_C_H2O(end,k) C_C_N2(end,k)];         % mole fraction of components in order: [C2H6 C2H4 O2 CO2 CO H2O N2]
    Initial_guess  = [C_gas_In C_Cpf_In C_Rof_In T0 Ts0];
    X              = fsolve(@(x) BoundaryEquations(x,Az,Ar,C_T(2:end-1,k),T0,u0,keffz,keffr,hw,Tb,epsilon,kg,hg,as,Density_bed,Flowin,Pt,R,RxnKinetic,Components,C_gas_In,C_gas_Out,'First','Length'),Initial_guess,Option);
    
    C_Cs_C2H6(1,k) = X(1);
    C_Cs_C2H4(1,k) = X(2);
    C_Cs_O2(1,k)   = X(3);
    C_Cs_CO2(1,k)  = X(4);
    C_Cs_CO(1,k)   = X(5);
    C_Cs_H2O(1,k)  = X(6);
    C_Cpf(1,k)     = X(7);
    C_Rof(1,k)     = X(8);
    C_T(1,k)       = X(9);
    C_Ts(1,k)      = X(10);
    
    %z=L
    C_C_C2H6(end,k)  = (((Az(end,2:end-1)*C_C_C2H6(2:end-1,k))+(Az(end,1)*C_C_C2H6(1,k)))/(-Az(end,end)));
    C_C_C2H4(end,k)  = (((Az(end,2:end-1)*C_C_C2H4(2:end-1,k))+(Az(end,1)*C_C_C2H4(1,k)))/(-Az(end,end)));
    C_C_O2(end,k)    = (((Az(end,2:end-1)*C_C_O2(2:end-1,k))+(Az(end,1)*C_C_O2(1,k)))/(-Az(end,end)));
    C_C_CO2(end,k)   = (((Az(end,2:end-1)*C_C_CO2(2:end-1,k))+(Az(end,1)*C_C_CO2(1,k)))/(-Az(end,end)));
    C_C_CO(end,k)    = (((Az(end,2:end-1)*C_C_CO(2:end-1,k))+(Az(end,1)*C_C_CO(1,k)))/(-Az(end,end)));
    C_C_H2O(end,k)   = (((Az(end,2:end-1)*C_C_H2O(2:end-1,k))+(Az(end,1)*C_C_H2O(1,k)))/(-Az(end,end)));
    C_C_N2(end,k)    = (((Az(end,2:end-1)*C_C_N2(2:end-1,k))+(Az(end,1)*C_C_N2(1,k)))/(-Az(end,end)));
    C_gas_In         = [C_C_C2H6(1,k) C_C_C2H4(1,k) C_C_O2(1,k) C_C_CO2(1,k) C_C_CO(1,k) C_C_H2O(1,k) C_C_N2(1,k)];         % mole fraction of components in order: [C2H6 C2H4 O2 CO2 CO H2O]
    C_gas_Out        = [C_C_C2H6(end-1,k) C_C_C2H4(end,k) C_C_O2(end,k) C_C_CO2(end,k) C_C_CO(end,k) C_C_H2O(end,k) C_C_N2(end,k)];         % mole fraction of components in order: [C2H6 C2H4 O2 CO2 CO H2O]
    Initial_guess    = [C_gas_Out C_Cpf_Out C_Rof_Out T0 Ts0];
    X                = fsolve(@(x) BoundaryEquations(x,Az,Ar,C_T(2:end-1,k),T0,u0,keffz,keffr,hw,Tb,epsilon,kg,hg,as,Density_bed,Flowin,Pt,R,RxnKinetic,Components,C_gas_In,C_gas_Out,'Last','Length'),Initial_guess,Option);
    
    C_Cs_C2H6(end,k) = X(1);
    C_Cs_C2H4(end,k) = X(2);
    C_Cs_O2(end,k)   = X(3);
    C_Cs_CO2(end,k)  = X(4);
    C_Cs_CO(end,k)   = X(5);
    C_Cs_H2O(end,k)  = X(6);
    C_Cpf(end,k)     = X(7);
    C_Rof(end,k)     = X(8);
    C_T(end,k)       = X(9);
    C_Ts(end,k)      = X(10);
    
end

for i = 1:numel(z_nodes)
    %r=0
    C_C_C2H6(i,1)  = (((Ar(1,2:end-1)*C_C_C2H6(i,2:end-1))-((Ar(1,end)/Ar(end,end))*Ar(end,2:end-1)*C_C_C2H6(i,2:end-1)))/(((Ar(1,end)/Ar(end,end))*Ar(end,1))-(Ar(1,1))));
    C_C_C2H4(i,1)  = (((Ar(1,2:end-1)*C_C_C2H4(i,2:end-1))-((Ar(1,end)/Ar(end,end))*Ar(end,2:end-1)*C_C_C2H4(i,2:end-1)))/(((Ar(1,end)/Ar(end,end))*Ar(end,1))-(Ar(1,1))));
    C_C_O2(i,1)    = (((Ar(1,2:end-1)*C_C_O2(i,2:end-1))-((Ar(1,end)/Ar(end,end))*Ar(end,2:end-1)*C_C_O2(i,2:end-1)))/(((Ar(1,end)/Ar(end,end))*Ar(end,1))-(Ar(1,1))));
    C_C_CO2(i,1)   = (((Ar(1,2:end-1)*C_C_CO2(i,2:end-1))-((Ar(1,end)/Ar(end,end))*Ar(end,2:end-1)*C_C_CO2(i,2:end-1)))/(((Ar(1,end)/Ar(end,end))*Ar(end,1))-(Ar(1,1))));
    C_C_CO(i,1)    = (((Ar(1,2:end-1)*C_C_CO(i,2:end-1))-((Ar(1,end)/Ar(end,end))*Ar(end,2:end-1)*C_C_CO(i,2:end-1)))/(((Ar(1,end)/Ar(end,end))*Ar(end,1))-(Ar(1,1))));
    C_C_H2O(i,1)   = (((Ar(1,2:end-1)*C_C_H2O(i,2:end-1))-((Ar(1,end)/Ar(end,end))*Ar(end,2:end-1)*C_C_H2O(i,2:end-1)))/(((Ar(1,end)/Ar(end,end))*Ar(end,1))-(Ar(1,1))));
    C_C_N2(i,1)    = (((Ar(1,2:end-1)*C_C_N2(i,2:end-1))-((Ar(1,end)/Ar(end,end))*Ar(end,2:end-1)*C_C_N2(i,2:end-1)))/(((Ar(1,end)/Ar(end,end))*Ar(end,1))-(Ar(1,1))));
    C_gas_In       = [C_C_C2H6(i,1) C_C_C2H4(i,1) C_C_O2(i,1) C_C_CO2(i,1) C_C_CO(i,1) C_C_H2O(i,1) C_C_N2(i,1)];         % mole fraction of components in order: [C2H6 C2H4 O2 CO2 CO H2O]
    C_gas_Out      = [C_C_C2H6(i,end) C_C_C2H4(i,end) C_C_O2(i,end) C_C_CO2(i,end) C_C_CO(i,end) C_C_H2O(i,end) C_C_N2(i,end)];         % mole fraction of components in order: [C2H6 C2H4 O2 CO2 CO H2O]
    Initial_guess  = [C_gas_In C_Cpf_In C_Rof_In T0 Ts0];
    X              = fsolve(@(x) BoundaryEquations(x,Az,Ar,C_T(i,2:end-1),T0,u0,keffz,keffr,hw,Tb,epsilon,kg,hg,as,Density_bed,Flowin,Pt,R,RxnKinetic,Components,C_gas_In,C_gas_Out,'First','Radius'),Initial_guess,Option);
    
    C_Cs_C2H6(i,1) = X(1);
    C_Cs_C2H4(i,1) = X(2);
    C_Cs_O2(i,1)   = X(3);
    C_Cs_CO2(i,1)  = X(4);
    C_Cs_CO(i,1)   = X(5);
    C_Cs_H2O(i,1)  = X(6);
    C_Cpf(i,1)     = X(7);
    C_Rof(i,1)     = X(8);
    C_T(i,1)       = X(9);
    C_Ts(i,1)      = X(10);
    
    %r=R
    C_C_C2H6(i,end)  = (((Ar(end,2:end-1)*C_C_C2H6(i,2:end-1))+(Ar(end,1)*(((Ar(1,2:end-1)*C_C_C2H6(i,2:end-1))-((Ar(1,end)/Ar(end,end))*Ar(end,2:end-1)*C_C_C2H6(i,2:end-1)))/(((Ar(1,end)/Ar(end,end))*Ar(end,1))-(Ar(1,1))))))/(-Ar(end,end)));
    C_C_C2H4(i,end)  = (((Ar(end,2:end-1)*C_C_C2H4(i,2:end-1))+(Ar(end,1)*(((Ar(1,2:end-1)*C_C_C2H4(i,2:end-1))-((Ar(1,end)/Ar(end,end))*Ar(end,2:end-1)*C_C_C2H4(i,2:end-1)))/(((Ar(1,end)/Ar(end,end))*Ar(end,1))-(Ar(1,1))))))/(-Ar(end,end)));
    C_C_O2(i,end)    = (((Ar(end,2:end-1)*C_C_O2(i,2:end-1))+(Ar(end,1)*(((Ar(1,2:end-1)*C_C_O2(i,2:end-1))-((Ar(1,end)/Ar(end,end))*Ar(end,2:end-1)*C_C_O2(i,2:end-1)))/(((Ar(1,end)/Ar(end,end))*Ar(end,1))-(Ar(1,1))))))/(-Ar(end,end)));
    C_C_CO2(i,end)   = (((Ar(end,2:end-1)*C_C_CO2(i,2:end-1))+(Ar(end,1)*(((Ar(1,2:end-1)*C_C_CO2(i,2:end-1))-((Ar(1,end)/Ar(end,end))*Ar(end,2:end-1)*C_C_CO2(i,2:end-1)))/(((Ar(1,end)/Ar(end,end))*Ar(end,1))-(Ar(1,1))))))/(-Ar(end,end)));
    C_C_CO(i,end)    = (((Ar(end,2:end-1)*C_C_CO(i,2:end-1))+(Ar(end,1)*(((Ar(1,2:end-1)*C_C_CO(i,2:end-1))-((Ar(1,end)/Ar(end,end))*Ar(end,2:end-1)*C_C_CO(i,2:end-1)))/(((Ar(1,end)/Ar(end,end))*Ar(end,1))-(Ar(1,1))))))/(-Ar(end,end)));
    C_C_H2O(i,end)   = (((Ar(end,2:end-1)*C_C_H2O(i,2:end-1))+(Ar(end,1)*(((Ar(1,2:end-1)*C_C_H2O(i,2:end-1))-((Ar(1,end)/Ar(end,end))*Ar(end,2:end-1)*C_C_H2O(i,2:end-1)))/(((Ar(1,end)/Ar(end,end))*Ar(end,1))-(Ar(1,1))))))/(-Ar(end,end)));
    C_C_N2(i,end)    = (((Ar(end,2:end-1)*C_C_N2(i,2:end-1))+(Ar(end,1)*(((Ar(1,2:end-1)*C_C_N2(i,2:end-1))-((Ar(1,end)/Ar(end,end))*Ar(end,2:end-1)*C_C_N2(i,2:end-1)))/(((Ar(1,end)/Ar(end,end))*Ar(end,1))-(Ar(1,1))))))/(-Ar(end,end)));
    C_gas_In         = [C_C_C2H6(i,1) C_C_C2H4(i,1) C_C_O2(i,1) C_C_CO2(i,1) C_C_CO(i,1) C_C_H2O(i,1) C_C_N2(i,1)];         % mole fraction of components in order: [C2H6 C2H4 O2 CO2 CO H2O]
    C_gas_Out        = [C_C_C2H6(i,end) C_C_C2H4(i,end) C_C_O2(i,end) C_C_CO2(i,end) C_C_CO(i,end) C_C_H2O(i,end) C_C_N2(i,end)];         % mole fraction of components in order: [C2H6 C2H4 O2 CO2 CO H2O]
    Initial_guess    = [C_gas_Out C_Cpf_Out C_Rof_Out T0 Ts0];
    X                = fsolve(@(x) BoundaryEquations(x,Az,Ar,C_T(i,2:end-1),T0,u0,keffz,keffr,hw,Tb,epsilon,kg,hg,as,Density_bed,Flowin,Pt,R,RxnKinetic,Components,C_gas_In,C_gas_Out,'Last','Radius'),Initial_guess,Option);
    
    C_Cs_C2H6(i,end) = X(1);
    C_Cs_C2H4(i,end) = X(2);
    C_Cs_O2(i,end)   = X(3);
    C_Cs_CO2(i,end)  = X(4);
    C_Cs_CO(i,end)   = X(5);
    C_Cs_H2O(i,end)  = X(6);
    C_Cpf(i,end)     = X(7);
    C_Rof(i,end)     = X(8);
    C_T(i,end)       = X(9);
    C_Ts(i,end)      = X(10);
end

for i = 1 : numel(z_nodes)
    for j = 1 : numel(r_nodes)
        C_C2H6_hat = C_C_C2H6(i,j) * Lz(i) * Lr(j) + C_C2H6_hat ;
        C_C2H4_hat = C_C_C2H4(i,j) * Lz(i) * Lr(j) + C_C2H4_hat ;
        C_O2_hat   = C_C_O2(i,j)   * Lz(i) * Lr(j) + C_O2_hat   ;
        C_CO2_hat  = C_C_CO2(i,j)  * Lz(i) * Lr(j) + C_CO2_hat  ;
        C_CO_hat   = C_C_CO(i,j)   * Lz(i) * Lr(j) + C_CO_hat   ;
        C_H2O_hat  = C_C_H2O(i,j)  * Lz(i) * Lr(j) + C_H2O_hat  ;
        C_N2_hat   = C_C_N2(i,j)   * Lz(i) * Lr(j) + C_N2_hat   ;
        Rof_hat    = C_Rof(i,j)    * Lz(i) * Lr(j) + Rof_hat    ;
        Cs_C2H6_hat = C_Cs_C2H6(i,j) * Lz(i) * Lr(j) + Cs_C2H6_hat ;
        Cs_C2H4_hat = C_Cs_C2H4(i,j) * Lz(i) * Lr(j) + Cs_C2H4_hat ;
        Cs_O2_hat   = C_Cs_O2(i,j)   * Lz(i) * Lr(j) + Cs_O2_hat   ;
        Cs_CO2_hat  = C_Cs_CO2(i,j)  * Lz(i) * Lr(j) + Cs_CO2_hat  ;
        Cs_CO_hat   = C_Cs_CO(i,j)   * Lz(i) * Lr(j) + Cs_CO_hat   ;
        Cs_H2O_hat  = C_Cs_H2O(i,j)  * Lz(i) * Lr(j) + Cs_H2O_hat  ;
        Cpf_hat     = C_Cpf(i,j)     * Lz(i) * Lr(j) + Cpf_hat     ;
        T_hat       = C_T(i,j)       * Lz(i) * Lr(j) + T_hat       ;
        Ts_hat      = C_Ts(i,j)      * Lz(i) * Lr(j) + Ts_hat      ;
    end
end
%% PLot

% Discretize the z and r direction
zz = linspace(0,1,20) ;
rr = linspace(0,1,10) ;

% Initialize the gas and solid concentration and temperature
C_C_C2H6  = zeros(numel(zz),numel(rr))    ;   C_C_C2H4  = zeros(numel(zz),numel(rr))  ;
C_C_O2    = zeros(numel(zz),numel(rr))    ;   C_C_CO2   = zeros(numel(zz),numel(rr))  ;
C_C_CO    = zeros(numel(zz),numel(rr))    ;   C_C_H2O   = zeros(numel(zz),numel(rr))  ;
C_C_N2    = zeros(numel(zz),numel(rr))    ;   C_Rof     = zeros(numel(zz),numel(rr))  ;
C_Cs_C2H6 = zeros(numel(zz),numel(rr))    ;   C_Cs_C2H4 = zeros(numel(zz),numel(rr))  ;
C_Cs_O2   = zeros(numel(zz),numel(rr))    ;   C_Cs_CO2  = zeros(numel(zz),numel(rr))  ;
C_Cs_CO   = zeros(numel(zz),numel(rr))    ;   C_Cs_H2O  = zeros(numel(zz),numel(rr))  ;
C_Cpf     = zeros(numel(zz),numel(rr))    ;   C_T       = zeros(numel(zz),numel(rr))  ;
C_Ts      = zeros(numel(zz),numel(rr))    ;

% Calculate the concentration and temperature for specified nodes in reactor
for i = 1 : numel(zz)
    for j = 1 : numel(rr)     
        C_C_C2H6(i,j)  = double(subs(C_C2H6_hat, [z,r], [zz(i),rr(j)]))    ;
        C_C_C2H4(i,j)  = double(subs(C_C2H4_hat, [z,r], [zz(i),rr(j)]))    ;
        C_C_O2(i,j)    = double(subs(C_O2_hat, [z,r], [zz(i),rr(j)]))      ;  
        C_C_CO2(i,j)   = double(subs(C_CO2_hat, [z,r], [zz(i),rr(j)]))     ;
        C_C_CO(i,j)    = double(subs(C_CO_hat, [z,r], [zz(i),rr(j)]))      ;
        C_C_H2O(i,j)   = double(subs(C_H2O_hat, [z,r], [zz(i),rr(j)]))     ;
        C_C_N2(i,j)    = double(subs(C_N2_hat, [z,r], [zz(i),rr(j)]))      ;
        C_Rof(i,j)     = double(subs(Rof_hat, [z,r], [zz(i),rr(j)]))       ;
        C_Cs_C2H6(i,j) = double(subs(Cs_C2H6_hat, [z,r], [zz(i),rr(j)]))   ;
        C_Cs_C2H4(i,j) = double(subs(Cs_C2H4_hat, [z,r], [zz(i),rr(j)]))   ;
        C_Cs_O2(i,j)   = double(subs(Cs_O2_hat, [z,r], [zz(i),rr(j)]))     ;
        C_Cs_CO2(i,j)  = double(subs(Cs_CO2_hat, [z,r], [zz(i),rr(j)]))    ;
        C_Cs_CO(i,j)   = double(subs(Cs_CO_hat, [z,r], [zz(i),rr(j)]))     ;
        C_Cs_H2O(i,j)  = double(subs(Cs_H2O_hat, [z,r], [zz(i),rr(j)]))    ;
        C_Cpf(i,j)     = double(subs(Cpf_hat, [z,r], [zz(i),rr(j)]))       ;
        C_T(i,j)       = double(subs(T_hat, [z,r], [zz(i),rr(j)]))         ;
        C_Ts(i,j)      = double(subs(Ts_hat, [z,r], [zz(i),rr(j)]))        ;
    end  
end

% Convert the normalize r,z direction to normal direction
Zz = (zmax - zmin)*zz + zmin ;
Rr = (rmax - rmin)*rr + rmin ;

% Initial concentration and Temperature of reactor
Initial_C_C2H6        =  ones(length(zz),length(rr))*C0_C2H6   ;
Initial_C_C2H4        =  ones(length(zz),length(rr))*C0_C2H4   ;
Initial_C_O2          =  ones(length(zz),length(rr))*C0_O2     ;
Initial_C_CO2         =  ones(length(zz),length(rr))*C0_CO2    ;
Initial_C_CO          =  ones(length(zz),length(rr))*C0_CO     ;
Initial_C_H2O         =  ones(length(zz),length(rr))*C0_H2O    ;
Initial_C_N2          =  ones(length(zz),length(rr))*C0_N2     ;

% conversion type(1) of [C2H6 C2H4 O2 CO2 CO] from equ (34),(35)

% Conv_C2H6 = abs(C_C_C2H6 - Initial_C_C2H6)./(Initial_C_C2H6)   ;
% Conv_C2H4 = (C_C_C2H4 - Initial_C_C2H4)./(Initial_C_C2H6)      ;
% Conv_O2   = abs(C_C_O2 - Initial_C_O2)./(Initial_C_O2)         ;
% Conv_CO2  = (C_C_CO2 - Initial_C_CO2)./(Initial_C_C2H6)        ;
% Conv_CO   = (C_C_CO - Initial_C_CO)./(Initial_C_C2H6)          ;

% conversion type(2) of [C2H6 C2H4 O2 CO2 CO] from equ (34),(35)

% Conv_C2H6 = abs(C_C_C2H6)./(Initial_C_C2H6)   ;
% Conv_C2H4 = (C_C_C2H4 - Initial_C_C2H4)./(Initial_C_C2H6)      ;
% Conv_O2   = abs(C_C_O2)./(Initial_C_O2)         ;
% Conv_CO2  = (C_C_CO2 - Initial_C_CO2)./(Initial_C_C2H6)        ;
% Conv_CO   = (C_C_CO - Initial_C_CO)./(Initial_C_C2H6)          ;

% conversion type(3) of [C2H6 C2H4 O2 CO2 CO] from equ (34),(35)

Conv_C2H6 = abs(C_C_C2H6)./(Initial_C_C2H6)   ;
Conv_C2H4 = (C_C_C2H4)./(Initial_C_C2H6)      ;
Conv_O2   = abs(C_C_O2)./(Initial_C_O2)         ;
Conv_CO2  = (C_C_CO2)./(Initial_C_C2H6)        ;
Conv_CO   = (C_C_CO)./(Initial_C_C2H6)          ;

% plot the converion and temperature in reactor length (considering with 
% mean of r direction conversion and temperature)

figure(1)

plot(Zz, mean(Conv_C2H6,2), '--g', Zz, mean(Conv_C2H4,2), 'c', ...
     Zz, mean(Conv_O2,2), '-.r', Zz, mean(Conv_CO2,2), '--b', ...
     Zz, mean(Conv_CO,2), ':m', 'LineWidth',2)
title('ODH reactor modeling - Conv & Yield')
xlabel('Reactor length [m]')
ylabel('Conversion and Yield')
grid minor

figure(2)

plot(Zz, mean(C_T,2), 'g', 'LineWidth',2)
title('ODH reactor modeling - Temp')
xlabel('Reactor length [m]')
ylabel('Temperature [oC]')
grid minor


%% Save the results
% save('Results.mat','Conv_C2H6','Conv_C2H4','Conv_O2','Conv_CO2','Conv_CO','C_T','Zz','Rr')


