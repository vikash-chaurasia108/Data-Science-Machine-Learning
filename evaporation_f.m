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


P8= 1.3053;

T1= 308.15;
T0= 298; 
n_p= 0.8;
n_t= 0.8;
n_m= 0.96;
del_TL= 10;
TH= 450;
mdot= 0.1;
y= 0.2;

R=	8.31451;
C0=	2.046009;
C1=	22.231991;
C2=	-11.658491;
C3=	2.691665;
Tc=	456.831;
MW=	152.93;


Tsat_8= interp1(P_sat,T_sat,P8)+273.15;
hg8= interp1(P_sat,hg,P8);
sg8=interp1(P_sat,sg,P8);
hf8= interp1(P_sat,hf,P8);
sf8=interp1(P_sat,sf,P8);

h1= hf8;
s1= sf8;


count=1;
for P6= 10:1:30;
    
    P6m(count)= P6;
    Tsat_6(count)= interp1(P_sat,T_sat,P6)+273.15;
    hg6(count)= interp1(P_sat,hg,P6);
    sg6(count)=interp1(P_sat,sg,P6);
    hf6(count)= interp1(P_sat,hf,P6);
    sf6(count)=interp1(P_sat,sf,P6);



 
    P7(count)= 3.73634*((P6).^0.03343);


    
    Tsat_7(count)= interp1(P_sat,T_sat,P7(count))+273.15;
    hg7(count)= interp1(P_sat,hg,P7(count));
    sg7(count)=interp1(P_sat,sg,P7(count));
    hf7(count)= interp1(P_sat,hf,P7(count));
    sf7(count)=interp1(P_sat,sf,P7(count));
    
    h9(count)= hf7(count);
    s9(count)= sf7(count);
    
    c7=1;
for T7m= Tsat_7(count):0.1:(Tsat_7(count)+5)
    T7{count}= T7m;
    
    delta7_1{count}= abs(sg6-(sg7(count)+(R/MW)*(C0*log(T7{count}./Tsat_7(count)) + ((C1/Tc)*(T7{count} - Tsat_7(count))) + (C2/(2*Tc^2))*(T7{count}.^2 - Tsat_7(count).^2) + (C3/(3*Tc^3))*(T7{count}.^3 - Tsat_7(count).^3))));
    [l7,m7]= cellfun(@min,delta7_1);
    
    T7s(count)= T7{count}(m7(count));
    c7= c7+1;
end

    h7s(count)= hg7(count)+(R/MW)*(C0*(T7s(count)-Tsat_7(count)) + ((C1/(2*Tc))*(T7s(count).^2 - Tsat_7(count)^2)) + (C2/(3*Tc^2))*(T7s(count).^3 - Tsat_7(count)^3) + (C3/(4*Tc^3))*(T7s(count).^4 - Tsat_7(count)^4));
    
    
c8=1;
for T8m= T1:0.1:(T1+15)
    T8(c8)= T8m;
    delta_8{count}= abs(sg6(count)-(sg8+(R/MW)*(C0*log(T8./Tsat_8) + ((C1/Tc)*(T8 - Tsat_8)) + (C2/(2*Tc^2))*(T8.^2 - Tsat_8^2) + (C3/(3*Tc^3))*(T8.^3 - Tsat_8^3))));
    [l,m]= cellfun(@min,delta_8);
    
    T8s= T8(m);
    c8= c8+1;
end

    h8xs(count)= hg8+(R/MW)*(C0*(T8s(count)-Tsat_8) + ((C1/(2*Tc))*(T8s(count).^2 - Tsat_8^2)) + (C2/(3*Tc^2))*(T8s(count).^3 - Tsat_8^3) + (C3/(4*Tc^3))*(T8s(count).^4 - Tsat_8^4));

    h8(count)= n_t*(h8xs(count)-hg6(count))+ hg6(count);

    
    c4=1;
for T4m= Tsat_7(count):0.1:(Tsat_7(count)+5)

    T4{count}(c4)= T4m;
    delta_4{count}= abs(s9(count)-(interp1(T_sat,sf, T4{count}-273.15)-P6.*0.0001));
    [l4,m4]= cellfun(@min,delta_4);
    
    T4s(count)=T4{count}(m4(count));
    c4=c4+1;
end

    h4s(count)= interp1(T_sat,hf,T4s(count)-273.15) + P6.*0.01;
    
    
     c2=1;
for T2m= T1:0.1:(T1+5)
    T2{count}(c2)= T2m;
    delta_2{count}= abs(s1-(interp1(T_sat,sf, T2{count}-273.15)-P6.*0.0001));
    [l2,m2]= cellfun(@min,delta_2);
    
    T2s(count)=T2{count}(m2(count));
    c2=c2+1;
end

    h2s(count)= interp1(T_sat,hf,T2s(count)-273.15) + P6.*0.01;
    
    
    count= count+1;
    
    a= count-1;
end

    
    
    h2= h1 + (h2s-h1)/n_p;
    h7= n_t*(h7s - hg6)+ hg6;
    
    h8s= h7-((h7 - h8)/n_t);
    
    h4= h9 + (h4s-h9)/n_p;

    h3= (y*(h7 - h9)/(1-y))+h2;
    
   
    
    h5= h3*(1-y)+ h4*y;
    
    
     
     
    Qdot= mdot*(hg6-h5);
   
    Wpdot= mdot*(y*(h4s-h9)+(1-y)*(h2s-h1))./n_p;
    
    Wtdot= mdot*((hg6-h7s)+(1-y)*(h7-h8s)).*n_t;
    
    Wnet= (Wtdot*n_m)- Wpdot;
    
   n_th= ((Wtdot*n_m)- Wpdot)./Qdot;


 
   n_orc= find(n_th == max(n_th(:)));
   
   nth= n_th(n_orc);
   
  n_II= n_th./(1-(T0/TH));
  
  TL= T1 - del_TL;
  
  Idot= T0*mdot*(((h5 - hg6)/TH) + (1-y)*((h8 - h1)/TL));

   
   figure;
hold on
title('System efficiency and pressure(cfoh)');
xlabel('Turbine Inlet Pressure(Bar)'); 
ylabel('efficiency');
xlim([8 33])
ylim([0 0.3])
plot(P6m,n_th)


figure;
hold on
title('Second Law efficiency and pressure(cfoh)');
xlabel('Turbine Inlet Pressure(Bar)'); 
ylabel('Second Law efficiency');
xlim([8 33])
ylim([0 1])
plot(P6m,n_II)


figure;
hold on
title('Total Exergy Destruction Rate and pressure(cfoh)');
xlabel('Turbine Inlet Pressure(Bar)'); 
ylabel('Total Exergy Destruction Rate(KJ/s');
plot(P6m, Idot)



    
    
    

