function F = BoundaryEquations(x,Az,Ar,C_T,T0,u0,keffz,keffr,hw,Tb,epsilon,kg,hg,as,Density_bed,Flowin,Pt,R,RxnKinetic,Components,CC_gas_In,CC_gas_Out,BoundaryCond,Direction)

CC_Cs_C2H6 = x(1);
CC_Cs_C2H4 = x(2);
CC_Cs_O2   = x(3);
CC_Cs_CO2  = x(4);
CC_Cs_CO   = x(5);
CC_Cs_H2O  = x(6);
CC_Cpf     = x(7);
CC_Rof     = x(8);
CC_T       = x(9);
CC_Ts      = x(10);

if strcmp(Direction,'Length')== 1
    
    if strcmp(BoundaryCond,'First')== 1
        
        CC_solid = CC_gas_In;
        CC_C_C2H6 = CC_gas_In(1);
        CC_C_C2H4 = CC_gas_In(2);
        CC_C_O2   = CC_gas_In(3);
        CC_C_CO2  = CC_gas_In(4);
        CC_C_CO   = CC_gas_In(5);
        CC_C_H2O  = CC_gas_In(6);
        CC_C_N2   = CC_gas_In(7);
        y_gas = CC_gas_In/sum(CC_gas_In);
        
        E_Cs_C2H6 = (1-epsilon)*kg*as*(CC_C_C2H6 - CC_Cs_C2H6) + Density_bed*(ODHReactions(CC_Cpf,CC_solid,CC_Ts,R,Pt,Flowin,RxnKinetic,[Components.deltaS0],[Components.deltaH0],1,'Mass'));
        E_Cs_C2H4 = (1-epsilon)*kg*as*(CC_C_C2H4 - CC_Cs_C2H4) + Density_bed*(ODHReactions(CC_Cpf,CC_solid,CC_Ts,R,Pt,Flowin,RxnKinetic,[Components.deltaS0],[Components.deltaH0],2,'Mass'));
        E_Cs_O2   = (1-epsilon)*kg*as*(CC_C_O2 - CC_Cs_O2) + Density_bed*(ODHReactions(CC_Cpf,CC_solid,CC_Ts,R,Pt,Flowin,RxnKinetic,[Components.deltaS0],[Components.deltaH0],3,'Mass'));
        E_Cs_CO2  = (1-epsilon)*kg*as*(CC_C_CO2 - CC_Cs_CO2) + Density_bed*(ODHReactions(CC_Cpf,CC_solid,CC_Ts,R,Pt,Flowin,RxnKinetic,[Components.deltaS0],[Components.deltaH0],4,'Mass'));
        E_Cs_CO   = (1-epsilon)*kg*as*(CC_C_CO - CC_Cs_CO) + Density_bed*(ODHReactions(CC_Cpf,CC_solid,CC_Ts,R,Pt,Flowin,RxnKinetic,[Components.deltaS0],[Components.deltaH0],5,'Mass'));
        E_Cs_H2O  = (1-epsilon)*kg*as*(CC_C_H2O - CC_Cs_H2O) + Density_bed*(ODHReactions(CC_Cpf,CC_solid,CC_Ts,R,Pt,Flowin,RxnKinetic,[Components.deltaS0],[Components.deltaH0],6,'Mass'));
        E_Cpf     = CC_Cpf - (((Components(1).cp_R(1) + Components(1).cp_R(2)*CC_T + Components(1).cp_R(3)*CC_T^2 + Components(1).cp_R(4)*CC_T^(-2))*R)*(y_gas(1)))*(1000/Components(1).Mw) ...
            - (((Components(2).cp_R(1) + Components(2).cp_R(2)*CC_T + Components(2).cp_R(3)*CC_T^2 + Components(2).cp_R(4)*CC_T^(-2))*R)*(y_gas(2)))*(1000/Components(2).Mw)...
            - (((Components(3).cp_R(1) + Components(3).cp_R(2)*CC_T + Components(3).cp_R(3)*CC_T^2 + Components(3).cp_R(4)*CC_T^(-2))*R)*(y_gas(3)))*(1000/Components(3).Mw) ...
            - (((Components(4).cp_R(1) + Components(4).cp_R(2)*CC_T + Components(4).cp_R(3)*CC_T^2 + Components(4).cp_R(4)*CC_T^(-2))*R)*(y_gas(4)))*(1000/Components(4).Mw) ...
            - (((Components(5).cp_R(1) + Components(5).cp_R(2)*CC_T + Components(5).cp_R(3)*CC_T^2 + Components(5).cp_R(4)*CC_T^(-2))*R)*(y_gas(5)))*(1000/Components(5).Mw) ...
            - (((Components(6).cp_R(1) + Components(6).cp_R(2)*CC_T + Components(6).cp_R(3)*CC_T^2 + Components(6).cp_R(4)*CC_T^(-2))*R)*(y_gas(6)))*(1000/Components(6).Mw) ...
            - (((Components(7).cp_R(1) + Components(7).cp_R(2)*CC_T + Components(7).cp_R(3)*CC_T^2 + Components(7).cp_R(4)*CC_T^(-2))*R)*(y_gas(7)))*(1000/Components(7).Mw);
        R = 8.2057e-5;                         % [(atm.m3)/(mol*K)]
        E_Rof = CC_Rof - 0.001*(Pt*((CC_C_C2H6/(CC_C_C2H6+CC_C_C2H4+CC_C_O2+CC_C_CO2+CC_C_CO+CC_C_H2O+CC_C_N2))*Components(1).Mw + (CC_C_C2H4/(CC_C_C2H6+CC_C_C2H4+CC_C_O2+CC_C_CO2+CC_C_CO+CC_C_H2O+CC_C_N2))*Components(2).Mw + (CC_C_O2/(CC_C_C2H6+CC_C_C2H4+CC_C_O2+CC_C_CO2+CC_C_CO+CC_C_H2O+CC_C_N2))*Components(3).Mw + (CC_C_CO2/(CC_C_C2H6+CC_C_C2H4+CC_C_O2+CC_C_CO2+CC_C_CO+CC_C_H2O+CC_C_N2))*Components(4).Mw + (CC_C_CO/(CC_C_C2H6+CC_C_C2H4+CC_C_O2+CC_C_CO2+CC_C_CO+CC_C_H2O+CC_C_N2))*Components(5).Mw ...
            + (CC_C_H2O/(CC_C_C2H6+CC_C_C2H4+CC_C_O2+CC_C_CO2+CC_C_CO+CC_C_H2O+CC_C_N2))*Components(6).Mw + (CC_C_N2/(CC_C_C2H6+CC_C_C2H4+CC_C_O2+CC_C_CO2+CC_C_CO+CC_C_H2O+CC_C_N2))*Components(7).Mw))/(R*CC_T); % [kg/m^3] Density of fluid;
        R=8.314e-3;                            % [(kJ/mol*K)]
        E_T       = CC_T - (u0*CC_Rof*CC_Cpf*T0 + keffz*Az(1,2:end-1)*C_T + (keffz*Az(1,end)/Az(end,end))*(-Az(end,2:end-1)*C_T))/((1 + (keffz*Az(end,1)*Az(1,end))/(Az(end,end)*(u0*CC_Rof*CC_Cpf - keffz*Az(1,1))))*(u0*CC_Rof*CC_Cpf - keffz*Az(1,1)));
        E_Ts      = (1-epsilon)*hg*as*(CC_T - CC_Ts) + Density_bed*ODHReactions(CC_Cpf,CC_solid,CC_Ts,R,Pt,Flowin,RxnKinetic,[Components.deltaS0],[Components.deltaH0],1,'Energy');
        
    elseif strcmp(BoundaryCond,'Last')== 1
        
        CC_solid = CC_gas_Out;
        CC_C_C2H6 = CC_gas_Out(1);
        CC_C_C2H4 = CC_gas_Out(2);
        CC_C_O2   = CC_gas_Out(3);
        CC_C_CO2  = CC_gas_Out(4);
        CC_C_CO   = CC_gas_Out(5);
        CC_C_H2O  = CC_gas_Out(6);
        CC_C_N2   = CC_gas_Out(7);
        y_gas = CC_gas_Out/sum(CC_gas_Out);
        
        E_Cs_C2H6 = (1-epsilon)*kg*as*(CC_C_C2H6 - CC_Cs_C2H6) + Density_bed*(ODHReactions(CC_Cpf,CC_solid,CC_Ts,R,Pt,Flowin,RxnKinetic,[Components.deltaS0],[Components.deltaH0],1,'Mass'));
        E_Cs_C2H4 = (1-epsilon)*kg*as*(CC_C_C2H4 - CC_Cs_C2H4) + Density_bed*(ODHReactions(CC_Cpf,CC_solid,CC_Ts,R,Pt,Flowin,RxnKinetic,[Components.deltaS0],[Components.deltaH0],2,'Mass'));
        E_Cs_O2   = (1-epsilon)*kg*as*(CC_C_O2 - CC_Cs_O2) + Density_bed*(ODHReactions(CC_Cpf,CC_solid,CC_Ts,R,Pt,Flowin,RxnKinetic,[Components.deltaS0],[Components.deltaH0],3,'Mass'));
        E_Cs_CO2  = (1-epsilon)*kg*as*(CC_C_CO2 - CC_Cs_CO2) + Density_bed*(ODHReactions(CC_Cpf,CC_solid,CC_Ts,R,Pt,Flowin,RxnKinetic,[Components.deltaS0],[Components.deltaH0],4,'Mass'));
        E_Cs_CO   = (1-epsilon)*kg*as*(CC_C_CO - CC_Cs_CO) + Density_bed*(ODHReactions(CC_Cpf,CC_solid,CC_Ts,R,Pt,Flowin,RxnKinetic,[Components.deltaS0],[Components.deltaH0],5,'Mass'));
        E_Cs_H2O  = (1-epsilon)*kg*as*(CC_C_H2O - CC_Cs_H2O) + Density_bed*(ODHReactions(CC_Cpf,CC_solid,CC_Ts,R,Pt,Flowin,RxnKinetic,[Components.deltaS0],[Components.deltaH0],6,'Mass'));
        E_Cpf     = CC_Cpf - (((Components(1).cp_R(1) + Components(1).cp_R(2)*CC_T + Components(1).cp_R(3)*CC_T^2 + Components(1).cp_R(4)*CC_T^(-2))*R)*(y_gas(1)))*(1000/Components(1).Mw) ...
            - (((Components(2).cp_R(1) + Components(2).cp_R(2)*CC_T + Components(2).cp_R(3)*CC_T^2 + Components(2).cp_R(4)*CC_T^(-2))*R)*(y_gas(2)))*(1000/Components(2).Mw) ...
            - (((Components(3).cp_R(1) + Components(3).cp_R(2)*CC_T + Components(3).cp_R(3)*CC_T^2 + Components(3).cp_R(4)*CC_T^(-2))*R)*(y_gas(3)))*(1000/Components(3).Mw) ...
            - (((Components(4).cp_R(1) + Components(4).cp_R(2)*CC_T + Components(4).cp_R(3)*CC_T^2 + Components(4).cp_R(4)*CC_T^(-2))*R)*(y_gas(4)))*(1000/Components(4).Mw) ...
            - (((Components(5).cp_R(1) + Components(5).cp_R(2)*CC_T + Components(5).cp_R(3)*CC_T^2 + Components(5).cp_R(4)*CC_T^(-2))*R)*(y_gas(5)))*(1000/Components(5).Mw) ...
            - (((Components(6).cp_R(1) + Components(6).cp_R(2)*CC_T + Components(6).cp_R(3)*CC_T^2 + Components(6).cp_R(4)*CC_T^(-2))*R)*(y_gas(6)))*(1000/Components(6).Mw) ...
            - (((Components(7).cp_R(1) + Components(7).cp_R(2)*CC_T + Components(7).cp_R(3)*CC_T^2 + Components(7).cp_R(4)*CC_T^(-2))*R)*(y_gas(7)))*(1000/Components(7).Mw);
        R = 8.2057e-5;                         % [(atm.m3)/(mol*K)]
        E_Rof     = CC_Rof - 0.001*(Pt*((CC_C_C2H6/(CC_C_C2H6+CC_C_C2H4+CC_C_O2+CC_C_CO2+CC_C_CO+CC_C_H2O+CC_C_N2))*Components(1).Mw + (CC_C_C2H4/(CC_C_C2H6+CC_C_C2H4+CC_C_O2+CC_C_CO2+CC_C_CO+CC_C_H2O+CC_C_N2))*Components(2).Mw + (CC_C_O2/(CC_C_C2H6+CC_C_C2H4+CC_C_O2+CC_C_CO2+CC_C_CO+CC_C_H2O+CC_C_N2))*Components(3).Mw + (CC_C_CO2/(CC_C_C2H6+CC_C_C2H4+CC_C_O2+CC_C_CO2+CC_C_CO+CC_C_H2O+CC_C_N2))*Components(4).Mw + (CC_C_CO/(CC_C_C2H6+CC_C_C2H4+CC_C_O2+CC_C_CO2+CC_C_CO+CC_C_H2O+CC_C_N2))*Components(5).Mw ...
            + (CC_C_H2O/(CC_C_C2H6+CC_C_C2H4+CC_C_O2+CC_C_CO2+CC_C_CO+CC_C_H2O+CC_C_N2))*Components(6).Mw + (CC_C_N2/(CC_C_C2H6+CC_C_C2H4+CC_C_O2+CC_C_CO2+CC_C_CO+CC_C_H2O+CC_C_N2))*Components(7).Mw))/(R*CC_T); % [kg/m^3] Density of fluid;
        R=8.314e-3;
        E_T       = CC_T - (((-Az(end,2:end-1)*C_T) - (Az(end,1)*((u0*CC_Rof*CC_Cpf*T0 + keffz*Az(1,2:end-1)*C_T + (keffz*Az(1,end)/Az(end,end))*(-Az(end,2:end-1)*C_T))/((1 + (keffz*Az(end,1)*Az(1,end))/(Az(end,end)*(u0*CC_Rof*CC_Cpf - keffz*Az(1,1))))*(u0*CC_Rof*CC_Cpf - keffz*Az(1,1))))))/(Az(end,end)));
        E_Ts      = (1-epsilon)*hg*as*(CC_T - CC_Ts) + Density_bed*ODHReactions(CC_Cpf,CC_solid,CC_Ts,R,Pt,Flowin,RxnKinetic,[Components.deltaS0],[Components.deltaH0],1,'Energy');
        
    end
elseif strcmp(Direction,'Radius')== 1
    
    if strcmp(BoundaryCond,'First')== 1
        
        CC_solid = CC_gas_In;
        CC_C_C2H6 = CC_gas_In(1);
        CC_C_C2H4 = CC_gas_In(2);
        CC_C_O2   = CC_gas_In(3);
        CC_C_CO2  = CC_gas_In(4);
        CC_C_CO   = CC_gas_In(5);
        CC_C_H2O  = CC_gas_In(6);
        CC_C_N2   = CC_gas_In(7);
        y_gas = CC_gas_In/sum(CC_gas_In);
        
        E_Cs_C2H6 = (1-epsilon)*kg*as*(CC_C_C2H6 - CC_Cs_C2H6) + Density_bed*(ODHReactions(CC_Cpf,CC_solid,CC_Ts,R,Pt,Flowin,RxnKinetic,[Components.deltaS0],[Components.deltaH0],1,'Mass'));
        E_Cs_C2H4 = (1-epsilon)*kg*as*(CC_C_C2H4 - CC_Cs_C2H4) + Density_bed*(ODHReactions(CC_Cpf,CC_solid,CC_Ts,R,Pt,Flowin,RxnKinetic,[Components.deltaS0],[Components.deltaH0],2,'Mass'));
        E_Cs_O2   = (1-epsilon)*kg*as*(CC_C_O2 - CC_Cs_O2) + Density_bed*(ODHReactions(CC_Cpf,CC_solid,CC_Ts,R,Pt,Flowin,RxnKinetic,[Components.deltaS0],[Components.deltaH0],3,'Mass'));
        E_Cs_CO2  = (1-epsilon)*kg*as*(CC_C_CO2 - CC_Cs_CO2) + Density_bed*(ODHReactions(CC_Cpf,CC_solid,CC_Ts,R,Pt,Flowin,RxnKinetic,[Components.deltaS0],[Components.deltaH0],4,'Mass'));
        E_Cs_CO   = (1-epsilon)*kg*as*(CC_C_CO - CC_Cs_CO) + Density_bed*(ODHReactions(CC_Cpf,CC_solid,CC_Ts,R,Pt,Flowin,RxnKinetic,[Components.deltaS0],[Components.deltaH0],5,'Mass'));
        E_Cs_H2O  = (1-epsilon)*kg*as*(CC_C_H2O - CC_Cs_H2O) + Density_bed*(ODHReactions(CC_Cpf,CC_solid,CC_Ts,R,Pt,Flowin,RxnKinetic,[Components.deltaS0],[Components.deltaH0],6,'Mass'));
        E_Cpf     = CC_Cpf - (((Components(1).cp_R(1) + Components(1).cp_R(2)*CC_T + Components(1).cp_R(3)*CC_T^2 + Components(1).cp_R(4)*CC_T^(-2))*R)*(y_gas(1)))*(1000/Components(1).Mw) ...
            - (((Components(2).cp_R(1) + Components(2).cp_R(2)*CC_T + Components(2).cp_R(3)*CC_T^2 + Components(2).cp_R(4)*CC_T^(-2))*R)*(y_gas(2)))*(1000/Components(2).Mw)...
            - (((Components(3).cp_R(1) + Components(3).cp_R(2)*CC_T + Components(3).cp_R(3)*CC_T^2 + Components(3).cp_R(4)*CC_T^(-2))*R)*(y_gas(3)))*(1000/Components(3).Mw) ...
            - (((Components(4).cp_R(1) + Components(4).cp_R(2)*CC_T + Components(4).cp_R(3)*CC_T^2 + Components(4).cp_R(4)*CC_T^(-2))*R)*(y_gas(4)))*(1000/Components(4).Mw) ...
            - (((Components(5).cp_R(1) + Components(5).cp_R(2)*CC_T + Components(5).cp_R(3)*CC_T^2 + Components(5).cp_R(4)*CC_T^(-2))*R)*(y_gas(5)))*(1000/Components(5).Mw) ...
            - (((Components(6).cp_R(1) + Components(6).cp_R(2)*CC_T + Components(6).cp_R(3)*CC_T^2 + Components(6).cp_R(4)*CC_T^(-2))*R)*(y_gas(6)))*(1000/Components(6).Mw) ...
            - (((Components(7).cp_R(1) + Components(7).cp_R(2)*CC_T + Components(7).cp_R(3)*CC_T^2 + Components(7).cp_R(4)*CC_T^(-2))*R)*(y_gas(7)))*(1000/Components(7).Mw);
        R = 8.2057e-5;                         % [(atm.m3)/(mol*K)]
        E_Rof = CC_Rof - 0.001*(Pt*((CC_C_C2H6/(CC_C_C2H6+CC_C_C2H4+CC_C_O2+CC_C_CO2+CC_C_CO+CC_C_H2O+CC_C_N2))*Components(1).Mw + (CC_C_C2H4/(CC_C_C2H6+CC_C_C2H4+CC_C_O2+CC_C_CO2+CC_C_CO+CC_C_H2O+CC_C_N2))*Components(2).Mw + (CC_C_O2/(CC_C_C2H6+CC_C_C2H4+CC_C_O2+CC_C_CO2+CC_C_CO+CC_C_H2O+CC_C_N2))*Components(3).Mw + (CC_C_CO2/(CC_C_C2H6+CC_C_C2H4+CC_C_O2+CC_C_CO2+CC_C_CO+CC_C_H2O+CC_C_N2))*Components(4).Mw + (CC_C_CO/(CC_C_C2H6+CC_C_C2H4+CC_C_O2+CC_C_CO2+CC_C_CO+CC_C_H2O+CC_C_N2))*Components(5).Mw ...
            + (CC_C_H2O/(CC_C_C2H6+CC_C_C2H4+CC_C_O2+CC_C_CO2+CC_C_CO+CC_C_H2O+CC_C_N2))*Components(6).Mw + (CC_C_N2/(CC_C_C2H6+CC_C_C2H4+CC_C_O2+CC_C_CO2+CC_C_CO+CC_C_H2O+CC_C_N2))*Components(7).Mw))/(R*CC_T); % [kg/m^3] Density of fluid;
        R=8.314e-3;                            % [(kJ/mol*K)]
        E_T       = CC_T - (-Ar(1,2:end-1)*C_T - Ar(1,end)*(keffr*Ar(end,2:end-1)*C_T + hw*Tb)/(hw - keffr*Ar(end,end)))/((1 + (keffr*Ar(1,end)*Ar(end,1))/((hw - keffr*Ar(end,end))*Ar(1,1)))*Ar(1,1));
        E_Ts      = (1-epsilon)*hg*as*(CC_T - CC_Ts) + Density_bed*ODHReactions(CC_Cpf,CC_solid,CC_Ts,R,Pt,Flowin,RxnKinetic,[Components.deltaS0],[Components.deltaH0],1,'Energy');
        
    elseif strcmp(BoundaryCond,'Last')== 1
        
        CC_solid = CC_gas_Out;
        CC_C_C2H6 = CC_gas_Out(1);
        CC_C_C2H4 = CC_gas_Out(2);
        CC_C_O2   = CC_gas_Out(3);
        CC_C_CO2  = CC_gas_Out(4);
        CC_C_CO   = CC_gas_Out(5);
        CC_C_H2O  = CC_gas_Out(6);
        CC_C_N2   = CC_gas_Out(7);
        y_gas = CC_gas_Out/sum(CC_gas_Out);
        
        E_Cs_C2H6 = (1-epsilon)*kg*as*(CC_C_C2H6 - CC_Cs_C2H6) + Density_bed*(ODHReactions(CC_Cpf,CC_solid,CC_Ts,R,Pt,Flowin,RxnKinetic,[Components.deltaS0],[Components.deltaH0],1,'Mass'));
        E_Cs_C2H4 = (1-epsilon)*kg*as*(CC_C_C2H4 - CC_Cs_C2H4) + Density_bed*(ODHReactions(CC_Cpf,CC_solid,CC_Ts,R,Pt,Flowin,RxnKinetic,[Components.deltaS0],[Components.deltaH0],2,'Mass'));
        E_Cs_O2   = (1-epsilon)*kg*as*(CC_C_O2 - CC_Cs_O2) + Density_bed*(ODHReactions(CC_Cpf,CC_solid,CC_Ts,R,Pt,Flowin,RxnKinetic,[Components.deltaS0],[Components.deltaH0],3,'Mass'));
        E_Cs_CO2  = (1-epsilon)*kg*as*(CC_C_CO2 - CC_Cs_CO2) + Density_bed*(ODHReactions(CC_Cpf,CC_solid,CC_Ts,R,Pt,Flowin,RxnKinetic,[Components.deltaS0],[Components.deltaH0],4,'Mass'));
        E_Cs_CO   = (1-epsilon)*kg*as*(CC_C_CO - CC_Cs_CO) + Density_bed*(ODHReactions(CC_Cpf,CC_solid,CC_Ts,R,Pt,Flowin,RxnKinetic,[Components.deltaS0],[Components.deltaH0],5,'Mass'));
        E_Cs_H2O  = (1-epsilon)*kg*as*(CC_C_H2O - CC_Cs_H2O) + Density_bed*(ODHReactions(CC_Cpf,CC_solid,CC_Ts,R,Pt,Flowin,RxnKinetic,[Components.deltaS0],[Components.deltaH0],6,'Mass'));
        E_Cpf     = CC_Cpf - (((Components(1).cp_R(1) + Components(1).cp_R(2)*CC_T + Components(1).cp_R(3)*CC_T^2 + Components(1).cp_R(4)*CC_T^(-2))*R)*(y_gas(1)))*(1000/Components(1).Mw) ...
            - (((Components(2).cp_R(1) + Components(2).cp_R(2)*CC_T + Components(2).cp_R(3)*CC_T^2 + Components(2).cp_R(4)*CC_T^(-2))*R)*(y_gas(2)))*(1000/Components(2).Mw)...
            - (((Components(3).cp_R(1) + Components(3).cp_R(2)*CC_T + Components(3).cp_R(3)*CC_T^2 + Components(3).cp_R(4)*CC_T^(-2))*R)*(y_gas(3)))*(1000/Components(3).Mw) ...
            - (((Components(4).cp_R(1) + Components(4).cp_R(2)*CC_T + Components(4).cp_R(3)*CC_T^2 + Components(4).cp_R(4)*CC_T^(-2))*R)*(y_gas(4)))*(1000/Components(4).Mw) ...
            - (((Components(5).cp_R(1) + Components(5).cp_R(2)*CC_T + Components(5).cp_R(3)*CC_T^2 + Components(5).cp_R(4)*CC_T^(-2))*R)*(y_gas(5)))*(1000/Components(5).Mw) ...
            - (((Components(6).cp_R(1) + Components(6).cp_R(2)*CC_T + Components(6).cp_R(3)*CC_T^2 + Components(6).cp_R(4)*CC_T^(-2))*R)*(y_gas(6)))*(1000/Components(6).Mw) ...
            - (((Components(7).cp_R(1) + Components(7).cp_R(2)*CC_T + Components(7).cp_R(3)*CC_T^2 + Components(7).cp_R(4)*CC_T^(-2))*R)*(y_gas(7)))*(1000/Components(7).Mw);
        R = 8.2057e-5;                         % [(atm.m3)/(mol*K)]
        E_Rof = CC_Rof - 0.001*(Pt*((CC_C_C2H6/(CC_C_C2H6+CC_C_C2H4+CC_C_O2+CC_C_CO2+CC_C_CO+CC_C_H2O+CC_C_N2))*Components(1).Mw + (CC_C_C2H4/(CC_C_C2H6+CC_C_C2H4+CC_C_O2+CC_C_CO2+CC_C_CO+CC_C_H2O+CC_C_N2))*Components(2).Mw + (CC_C_O2/(CC_C_C2H6+CC_C_C2H4+CC_C_O2+CC_C_CO2+CC_C_CO+CC_C_H2O+CC_C_N2))*Components(3).Mw + (CC_C_CO2/(CC_C_C2H6+CC_C_C2H4+CC_C_O2+CC_C_CO2+CC_C_CO+CC_C_H2O+CC_C_N2))*Components(4).Mw + (CC_C_CO/(CC_C_C2H6+CC_C_C2H4+CC_C_O2+CC_C_CO2+CC_C_CO+CC_C_H2O+CC_C_N2))*Components(5).Mw ...
            + (CC_C_H2O/(CC_C_C2H6+CC_C_C2H4+CC_C_O2+CC_C_CO2+CC_C_CO+CC_C_H2O+CC_C_N2))*Components(6).Mw + (CC_C_N2/(CC_C_C2H6+CC_C_C2H4+CC_C_O2+CC_C_CO2+CC_C_CO+CC_C_H2O+CC_C_N2))*Components(7).Mw))/(R*CC_T); % [kg/m^3] Density of fluid;
        R=8.314e-3;                            % [(kJ/mol*K)]
        E_T       = CC_T - (keffr*Ar(end,2:end-1)*C_T + hw*Tb + keffr*Ar(end,1)*((-Ar(1,2:end-1)*C_T - Ar(1,end)*(keffr*Ar(end,2:end-1)*C_T + hw*Tb)/(hw - keffr*Ar(end,end)))/((1 + (keffr*Ar(1,end)*Ar(end,1))/((hw - keffr*Ar(end,end))*Ar(1,1)))*Ar(1,1))))/(hw - keffr*Ar(end,end));
        E_Ts      = (1-epsilon)*hg*as*(CC_T - CC_Ts) + Density_bed*ODHReactions(CC_Cpf,CC_solid,CC_Ts,R,Pt,Flowin,RxnKinetic,[Components.deltaS0],[Components.deltaH0],1,'Energy');
        
    end
end

F = [E_Cs_C2H6 E_Cs_C2H4 E_Cs_O2 E_Cs_CO2 E_Cs_CO E_Cs_H2O E_Cpf E_Rof E_T E_Ts];

end