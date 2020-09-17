function rxn = ODHReactions(Cpf,CC_s,Ts,R,Pt,Flowin,RxnKinetic,deltaS0,deltaH0,compnumber,type)
% This code is for modeling of ODH reaction kinetics
Tstar=400; %************** Assumption [=] C **************
% n_solid  = Flowin * CC_s ; % mole flow of each component [Nm^3/s * mol/m^3] = [mol/s]
% n_react = n_solid;         % mole of products in reaction
% nt_solid = sum(n_solid);     [unused]
Ct_solid = sum(CC_s);   % total mole concentration in solid phase

% component order list: [C2H6 C2H4 O2 CO2 CO H2O N2]
P_solid = Pt*(CC_s/Ct_solid);

% component order list for reaction: [C2H6 C2H4 O2 CO2 CO H2O]
K_O2 = exp((deltaS0(3) - deltaH0(3)*(1/Ts - 1/(Tstar+273.15)))/R);
K_H2O = exp((deltaS0(6) - deltaH0(6)*(1/Ts - 1/(Tstar+273.15)))/R);
Tetha_star = 1/(1 + (K_O2*P_solid(3))^0.5 + K_H2O*P_solid(6));
Tetha_O    = sqrt(K_O2*P_solid(3))*Tetha_star;

k = ones(1,5);
rxn = k(:);
for i = 1:5
    k(i)   = exp(RxnKinetic.Aprime(i) - RxnKinetic.EnergyA(i)/R*(1/Ts - 1/(Tstar+273.15)));
    if i < 4
        rxn(i) = k(i) * P_solid(1) * Tetha_O^(RxnKinetic.m(i));
    else
        rxn(i) = k(i) * P_solid(2) * Tetha_O^(RxnKinetic.m(i));
    end
end


if strcmp(type,'Energy') == 1
    %     n_product = n_solid;
    %     sum of heat reactions
    %     for i = 1:6
    %         n_product(i) = n_react(i) + rxn' * RxnKinetic.vcoffrxn(:,i); % * mcat must be consider
    %     end
    %     deltaH_reactants = sum(n_react * Cpf * (298.15 - Ts));
    %     deltaH_products  = sum(n_product * Cpf * (Ts - 298.15));
    deltaH_std = [RxnKinetic.deltaHstd]';
    rxn = sum(rxn' * (deltaH_std));% + 0*deltaH_reactants + 0*deltaH_products));
elseif strcmp(type,'Mass') == 1
    % sum of rate of reactions
    rxn = rxn' * RxnKinetic.vcoffrxn(:,compnumber);
end
end