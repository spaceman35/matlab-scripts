close all
clear all
clc
%% Initialization of variables

%%% Decleration of Variables with Initial Condition %%%%% We need to decleare these variables in float %%%%%%
r=2.0e7;   v_r=-6000;  v_theta= 5000;  theta=0;                            
d_t= 5;                             
r_vec = []; theta_vec = []; R_vec = []; Theta_vec = [];
a=2*r;   %%% supposing r to be rp 

%%% Initilaization of variable
G=6.67408e-11;                   
Mp=5.972e24 ;                    
dr_n = v_r *d_t;                     
dtheta_n= (v_theta * d_t)/r ;    
r_n=r;                           
theta_n=theta;                   

%%%%% Using Initial Conditions %%%%%%
  
%%%%% Caclutaion of Time required to stay in the desired elliptical orbit %%%%%%%%
mu = G*Mp; 
T_Assumed= (2*pi*sqrt(a*a*a))/(sqrt(mu)) %%%assumed Time period
 
 %%% Writing The recursive Equations%%%%%%
%% Cacluation of Recursive Values
for i=1:d_t:T_Assumed

    %%%%Calculating rdaial distance and angle of orbit
     r_n1= r_n + dr_n; 
    theta_n1= theta_n + dtheta_n;
 
    dr_n1 = dr_n+(((r_n+(0.5*dr_n))*(dtheta_n^2)) - (mu /(r_n^2)*(d_t^2)));
    dtheta_n1 = dtheta_n - [(2*dr_n*dtheta_n)/(r_n + (0.5*dr_n))];
  
    %%% Assigning the delta values %%%%%
   
   dr_n = dr_n1;    
   dtheta_n = dtheta_n1;
   r_n = r_n1;
   theta_n= theta_n1;
   
   %%%%% Creating the Distance Vector and Angle Vector
  r_vec = [r_vec r_n1];
  theta_vec = [theta_vec theta_n1];
end
%% Calculation of Ecentricity and Time period
e= (max(r_vec)-r)/(max(r_vec)+r)
sma = (max(r_vec)+r)/2;
T_Aprrox= (2*pi*sqrt(sma^3))/(sqrt(mu))


%% Plotting the calculated distance and Angle vectors in rectangular and polar cordinates 

figure;
   subplot (1,3,1); plot(r_vec);title('Distnace Vector with Vr= -6000');
   subplot (1,3,2); plot(theta_vec); title('Angle Vcetor with Vr= -6000');
   subplot (1,3,3); polarplot(theta_vec,r_vec); title('Orbital plot in Polar Cordinates with Vr= -6000');

 %% Recalculating the vectors based on Cacluated Time value
 %%%nitializing the values again 
DR_n = v_r *d_t;                     
DTheta_n= (v_theta * d_t)/r ;    
R_n=r;                           
Theta_n=theta; 

 for j=1:d_t:T_Aprrox

    %%%%Calculating rdaial distance and angle of orbit
     R_n1= R_n + DR_n; 
    Theta_n1= Theta_n + DTheta_n;
    
    DR_n1 = DR_n+(((R_n+(0.5*DR_n))*(DTheta_n^2)) - (mu /(R_n^2)*(d_t^2)));
    DTheta_n1 = DTheta_n - [(2*DR_n*DTheta_n)/(R_n + (0.5*DR_n))];
   
  
    %%% Assigning the delta values %%%%%
   
   DR_n = DR_n1;
   DTheta_n = DTheta_n1;
   R_n = R_n1;
   Theta_n= Theta_n1;
   
   %%%%% Creating the Distance Vector and Angle Vector
  R_vec = [R_vec R_n1];
  Theta_vec = [Theta_vec Theta_n1];
 end
 %%%%%% Cacluating The actual Time Period %%%%%%%
 
 
 %% Replotting the values
figure;
   subplot (1,3,1); plot(R_vec);title('Actual Distnace Vector with Vr= -6000');
   subplot (1,3,2); plot(Theta_vec); title('Actual Angle Vcetor with Vr= -6000');
   subplot (1,3,3); polarplot(Theta_vec,R_vec); title('Actual Orbital plot in Polar Cordinates with Vr= -6000');


   