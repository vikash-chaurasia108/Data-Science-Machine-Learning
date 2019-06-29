clear all
clc

[~, ~, raw] = xlsread('C:\Users\user\Desktop\sat.xlsx','R123','A2:F185');

data = reshape([raw{:}],size(raw));

sat = table;

sat.Tsat = data(:,1);
sat.P = data(:,2);
sat.hf = data(:,3);
sat.hg = data(:,4);
sat.sf = data(:,5);
sat.sg = data(:,6);

clearvars data raw;

format long

T_sat= sat{:,1};
P_sat= sat{:,2};
hf= sat{:,3};
hg= sat{:,4};
sf= sat{:,5};
sg= sat{:,6};

P6= 10;
P8= 1.3053;
P7= 4.0353;
T1= 308.15;
T0= 298; 
n_p= 0.8;
n_t= 0.8;
n_m= 0.96;
del_TL= 10;
TH= 425;
mdot= 0.1;
y= 0.2;

R=	8.31451;
C0=	2.046009;
C1=	22.231991;
C2=	-11.658491;
C3=	2.691665;
Tc=	456.831;
MW=	152.93;

Tsat_6= interp1(P_sat,T_sat,P6)+273.15;
hg6= interp1(P_sat,hg,P6);
sg6=interp1(P_sat,sg,P6);
hf6= interp1(P_sat,hf,P6);
sf6=interp1(P_sat,sf,P6);


Tsat_8= interp1(P_sat,T_sat,P8)+273.15;
hg8= interp1(P_sat,hg,P8);
sg8=interp1(P_sat,sg,P8);
hf8= interp1(P_sat,hf,P8);
sf8=interp1(P_sat,sf,P8);


Tsat_7= interp1(P_sat,T_sat,P7)+273.15;
hg7= interp1(P_sat,hg,P7);
sg7=interp1(P_sat,sg,P7);
hf7= interp1(P_sat,hf,P7);
sf7=interp1(P_sat,sf,P7);

count=1;
for T= Tsat_6:1:(TH-10)
    
    TIT(count)= T;
    
    h6(count)= hg6+(R/MW)*(C0*(T-Tsat_6) + ((C1/(2*Tc))*(T^2 - Tsat_6^2)) + (C2/(3*Tc^2))*(T^3 - Tsat_6^3) + (C3/(4*Tc^3))*(T^4 - Tsat_6^4));
    s6(count)= sg6+(R/MW)*(C0*log(T/Tsat_6) + ((C1/Tc)*(T - Tsat_6)) + (C2/(2*Tc^2))*(T^2 - Tsat_6^2) + (C3/(3*Tc^3))*(T^3 - Tsat_6^3));
    

c8=1;
for T8m= T1:0.1:(T1+200)
    T8(c8)= T8m;
    delta_8{count}= abs(s6(count)-(sg8+(R/MW)*(C0*log(T8./Tsat_8) + ((C1/Tc)*(T8 - Tsat_8)) + (C2/(2*Tc^2))*(T8.^2 - Tsat_8^2) + (C3/(3*Tc^3))*(T8.^3 - Tsat_8^3))));
    [l,m]= cellfun(@min,delta_8);
    
    T8s=T8(m);
    c8= c8+1;
end

    

    h8xs= hg8+(R/MW)*(C0*(T8s-Tsat_8) + ((C1/(2*Tc))*(T8s.^2 - Tsat_8^2)) + (C2/(3*Tc^2))*(T8s.^3 - Tsat_8^3) + (C3/(4*Tc^3))*(T8s.^4 - Tsat_8^4));

    
    
     
    c7=1;
for T7m= Tsat_7:0.1:(Tsat_7+200)
    T7(c7)= T7m;
    delta_7{count}= abs(s6(count)-(sg7+(R/MW)*(C0*log(T7./Tsat_7) + ((C1/Tc)*(T7 - Tsat_7)) + (C2/(2*Tc^2))*(T7.^2 - Tsat_7^2) + (C3/(3*Tc^3))*(T7.^3 - Tsat_7^3))));
    [l7,m7]= cellfun(@min,delta_7);
    
    T7s=T7(m7);
    c7= c7+1;
end
    

    h7s= hg7+(R/MW)*(C0*(T7s - Tsat_7) + ((C1/(2*Tc))*(T7s.^2 - Tsat_7^2)) + (C2/(3*Tc^2))*(T7s.^3 - Tsat_7^3) + (C3/(4*Tc^3))*(T7s.^4 - Tsat_7^4));
    
    
    count= count+1;
end

    h8=  h6- n_t*(h6-h8xs);

    h7= hg6- n_t*(hg6-h7s);
    
    h8s= h7-((h7-h8)./n_t);
    
    
    
    
    h1= hf8;
    s1= sf8;
      
    c2=1;
for T2m= T1:0.01:(T1+5)
    Tl(c2)= T2m;
    delta_2(c2)= abs(s1-(interp1(T_sat,sf,T2m-273.15)-P6*0.0001));
    g2= find(delta_2 == min(delta_2(:)));
    
    T2s=Tl(g2);
    c2= c2+1;
end

    h2s= interp1(T_sat,hf,T2s-273.15) + P6*0.01;
    h2= h1 + (h2s-h1)/n_p;
    
    
    
    
    h9= hf7;
    s9= sf7;
    
     c4=1;
for T4m= Tsat_7:0.01:(Tsat_7+5)
    T4(c4)= T4m;
    delta_4(c4)= abs(s9-(interp1(T_sat,sf,T4m-273.15)-P6*0.0001));
    g4= find(delta_4 == min(delta_4(:)));
    
    T4s=T4(g4);
    c4= c4+1;
end

    h4s= interp1(T_sat,hf,T4s-273.15) + P6*0.01;
    h4= h9 + (h4s-h9)/n_p;
    
    h3= (y*(h7-h9)/(1-y))+h2;
    
    h5= h3*(1-y)+ h4*y;
    
    
    Qdot= mdot*(h6-h5);
   
    Wpdot= mdot*(y*(h4s-h9)+(1-y)*(h2s-h1))./n_p;
    
    Wtdot= mdot*((h6-h7s)+(1-y)*(h7-h8s)).*n_t;
    
    Wnet= (Wtdot*n_m)- Wpdot;
    
    n_th= (n_m*Wtdot- Wpdot)./Qdot;
    
    n_II= n_th./(1-(T0./TH));
    
    TL= T1 - del_TL;
    
    Idot= T0*mdot*(((h5 - h6)/TH) + (1-y)*((h8 - h1)/TL));

    
    figure;
hold on
title('System efficiency and the availability ratio versus turbine inlet temperature');
xlabel('Turbine inlet temperature(K)'); 
ylabel('efficiency');
ylim([0 0.3])
plot(TIT,n_th)

figure;
hold on
title('Variation of the second law efficiency with turbine inlet temperature');
xlabel('Turbine inlet temperature(K)'); 
ylabel('Second law efficiency)');
ylim([0 0.7])
plot(TIT,n_II)
    
figure;
hold on
title('Total Exergy Destruction Rate and turbine inlet temperature');
xlabel('Turbine Inlet Temperature(k)'); 
ylabel('Total Exergy Destruction Rate(KJ/s');
plot(TIT, Idot)


