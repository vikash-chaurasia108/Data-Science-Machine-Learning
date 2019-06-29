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
T1= 308;
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

  c8=1;
for T8m= T1:0.01:T1+5
    T8(c8)= T8m;
    delta_8(c8)= abs(sg6-(sg8+(R/MW)*(C0*log(T8m/Tsat_8) + ((C1/Tc)*(T8m - Tsat_8)) + (C2/(2*Tc^2))*(T8m.^2 - Tsat_8^2) + (C3/(3*Tc^3))*(T8m.^3 - Tsat_8^3))));
    g8= find(delta_8 == min(delta_8(:)));
    
    T8s=T8(g8);
    c8= c8+1;
end

    h6= hg6;
    s6= sg6;

    h8xs= hg8+(R/MW)*(C0*(T8s-Tsat_8) + ((C1/(2*Tc))*(T8s.^2 - Tsat_8^2)) + (C2/(3*Tc^2))*(T8s.^3 - Tsat_8^3) + (C3/(4*Tc^3))*(T8s.^4 - Tsat_8^4));

    h8= n_t*(h8xs-h6)+ h6;
    
    
    h1= hf8;
    s1= sf8;
    
    
    
    c2=1;
for T2m= T1:0.01:T1+5
    Tl(c2)= T2m;
    delta_2(c2)= abs(s1-(interp1(T_sat,sf,T2m-273.15)-P6*0.0001));
    g2= find(delta_2 == min(delta_2(:)));
    
    T2s=Tl(g2);
    c2= c2+1;
end

    h2s= interp1(T_sat,hf,T2s-273.15) + P6*0.01;
    h2= h1 + (h2s-h1)/n_p;
    
    count=1
 for P7x= P8:0.01:P6
    
    P_cfoh(count)=P7x; 
    Tsat_7x(count)= interp1(P_sat,T_sat,P7x)+273.15;
    hg7x(count)= interp1(P_sat,hg,P7x);
    sg7x(count)=interp1(P_sat,sg,P7x);
    hf7x(count)= interp1(P_sat,hf,P7x);
    sf7x(count)=interp1(P_sat,sf,P7x);
    
    
   c7=1;
for T7m= Tsat_7x(count):0.1:Tsat_7x(count)+5
    T7{count}(c7)= T7m;
    
    delta7_1{count}= abs(sg6-(sg7x(count)+(R/MW)*(C0*log(T7{count}/Tsat_7x(count)) + ((C1/Tc)*(T7{count} - Tsat_7x(count))) + (C2/(2*Tc^2))*(T7{count}.^2 - Tsat_7x(count).^2) + (C3/(3*Tc^3))*(T7{count}.^3 - Tsat_7x(count).^3))));
    [l7,m7]= cellfun(@min,delta7_1);
    
    T7xs(count)=T7{count}(m7(count));
    c7= c7+1;
end
    h7xs(count)= hg7x(count)+(R/MW)*(C0*(T7xs(count)-Tsat_7x(count)) + ((C1/(2*Tc))*(T7xs(count).^2 - Tsat_7x(count)^2)) + (C2/(3*Tc^2))*(T7xs(count).^3 - Tsat_7x(count)^3) + (C3/(4*Tc^3))*(T7xs(count).^4 - Tsat_7x(count)^4));
    
   
    
    
    
    h9(count)= hf7x(count);
    s9(count)= sf7x(count);
    
    
    
    
    c4=1;
for T4m= Tsat_7x(count):0.1:Tsat_7x(count)+5
    T4{count}(c4)= T4m;
    delta_4{count}= abs(s9(count)-(interp1(T_sat,sf, T4{count}-273.15)-P6*0.0001));
    [l4,m4]= cellfun(@min,delta_4);
    
    T4s(count)=T4{count}(m4(count));
    c4=c4+1;
end

    h4s(count)= interp1(T_sat,hf,T4s(count)-273.15) + P6*0.01;
    
   

   count=count+1;
   a=count-1;
 end
 
    
    
     h7x= n_t*(h7xs-h6)+ h6;
    
    h8s= h7x-((h7x-h8)/n_t);
    
    h4= h9 + (h4s-h9)/n_p;

    h3= (y*(h7x-h9)/(1-y))+h2;
    
   
    
    h5= h3*(1-y)+ h4*y;
    
    c=1;
  for count= 1:1:a 
    if (h3(count) <= h9(count)&& h7x(count) <= h6)
     
    Qdot= mdot*(h6-h5);
   
    Wpdot= mdot*(y*(h4s-h9)+(1-y)*(h2s-h1))./n_p;
    
    Wtdot= mdot*((h6-h7xs)+(1-y)*(h7x-h8s)).*n_t;
    
    
    
   n_th(count)= ((Wtdot(count)*n_m)- Wpdot(count))/Qdot(count);
   
  
   
   

    else n_th(count) = 0;

  end
 c= c+1;
  end
  
   n_orc= find(n_th == max(n_th(:)));
   
   nth= n_th(n_orc);
  
  P7= P_cfoh(n_orc);
   
   
   figure;
hold on
title('System efficiency and pressure(cfoh)');
xlabel('pressure(cfoh)(Bar)'); 
ylabel('efficiency');
ylim([0 0.3])
plot(P_cfoh,n_th)