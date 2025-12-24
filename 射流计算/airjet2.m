%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%% Properties of the jet %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The properties of the jet are entered here
clear all;
close all;
clc;

Ujet = 1;    
R = 0.01;     
H = 60*R;     
lr = H/R;    
b = H-R;  

delta_G = 0.5*R;  

Dr = 769;	

Vr = 55.6;

qqq=10;  %气相速度
Uww=1;  %液相速度
kkk=3;  %K

Ugstar = 1.2*Uww;	

Reee = qqq*R*1.25/1.79*100000;
Weee = 1.25*qqq^2*R/0.0072;
Rel_matrix = [Reee; Reee];	% Reynolds number using liquid viscosity, Ujet and radius of the jet (There are 2 values which are 
                   
Wel_matrix = [Weee; Weee];    % Weber number using liquid density, Ujet, radius of the jet and the surface tension coefficient, gamma
                   


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%% Prompt for the number of divisions (cells) %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% In this section, number of divisions (cells) in liquid and gas domains
% will be entered separately. Default number of divisions in liquid (Nl=40)
% and default number of divisions in gas (Ng=70) will be used. The code can
% be run for different values of Nl and Ng. It should be made sure that
% sufficient number of divisions are used such that the dispersion curves
% converge as Nl and Ng are increased.

%% Number of divisions (cells) in liquid

Nl = 30;

%% Number of divisions (cells) in gas
Ng = 60;


%% α输入区间
ast =  0.1;        
mmm = 1/R/100; 
aend = 1000*mmm;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%% Defining Chebyshev matrices in different domains %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Liquid part
% Chebyshev matrices are defined for liquid domain except for
% liquid centerline (r=0 or \tilde{r}_L=-1} and
% interface on liquid side (r=R or \tilde{r}_L=+1)

for j=1:Nl-1
    rL(j,1) = -cos((pi/((Nl)))*j);    % rL=(-1,+1).
end

vec_liq = acos(rL);

% The Chebyshev matrices are defined below 
D0_liq=[];
vec=(1:1:Nl-1)';
for j=0:1:Nl
    D0_liq=[D0_liq cos(j*vec_liq)];
end
lv_liq=length(vec);
D1_liq=[zeros(lv_liq,1) D0_liq(:,1) 4*D0_liq(:,2)];
D2_liq=[zeros(lv_liq,1) zeros(lv_liq,1) 4*D0_liq(:,1)];
for j=3:Nl
    D1_liq=[D1_liq 2*j*D0_liq(:,j)+j*D1_liq(:,j-1)/(j-2)];
    D2_liq=[D2_liq 2*j*D1_liq(:,j)+j*D2_liq(:,j-1)/(j-2)];
end

%% Gas part
% Chebyshev matrices are defined for gas domain except for
% gas boundary (r=L or \tilde{r}_G=+1} and
% interface on gas side (r=R or \tilde{r}_G=-1)

for j=1:Ng-1
    rG(j,1) = -cos((pi/((Ng)))*j);  % rG=(-1,+1).
end
vec_gas = acos(rG);

% The Chebyshev matrices are defined below 
D0_gas=[];
for j=0:1:Ng
    D0_gas=[D0_gas cos(j*vec_gas)];
end
lv_gas=length(vec_gas);
D1_gas=[zeros(lv_gas,1) D0_gas(:,1) 4*D0_gas(:,2)];
D2_gas=[zeros(lv_gas,1) zeros(lv_gas,1) 4*D0_gas(:,1)];
for j=3:Ng
    D1_gas=[D1_gas 2*j*D0_gas(:,j)+j*D1_gas(:,j-1)/(j-2)];
    D2_gas=[D2_gas 2*j*D1_gas(:,j)+j*D2_gas(:,j-1)/(j-2)];
end

%% Liquid centerline
% Chebyshev matrices are defined for liquid for
% liquid centerline (r=0 or \tilde{r}_L=-1}

vec_liqCenter = pi;   
rL_liqCenter = cos(vec_liqCenter); % rL=-1.
    
% The Chebyshev matrices are defined below
D0_liqCenter=[];
for j=0:1:Nl
    D0_liqCenter=[D0_liqCenter cos(j*vec_liqCenter)];
end
lv_liqCenter=length(vec_liqCenter);
D1_liqCenter=[zeros(lv_liqCenter,1) D0_liqCenter(:,1) 4*D0_liqCenter(:,2)];
D2_liqCenter=[zeros(lv_liqCenter,1) zeros(lv_liqCenter,1) 4*D0_liqCenter(:,1)];
for j=3:Nl
    D1_liqCenter=[D1_liqCenter 2*j*D0_liqCenter(:,j)+j*D1_liqCenter(:,j-1)/(j-2)];
    D2_liqCenter=[D2_liqCenter 2*j*D1_liqCenter(:,j)+j*D2_liqCenter(:,j-1)/(j-2)];
end
    
%% Gas boundary
% Chebyshev matrices are defined for 
% gas boundary (r=L or \tilde{r}_G=+1}

vec_gasBoundary = pi;   
rG_gasBoundary = cos(vec_gasBoundary);% rG=-1.
    
% The Chebyshev matrices are defined below
D0_gasBoundary=[];
for j=0:1:Ng
    D0_gasBoundary=[D0_gasBoundary cos(j*vec_gasBoundary)];
end
lv_gasBoundary=length(vec_gasBoundary);
D1_gasBoundary=[zeros(lv_gasBoundary,1) D0_gasBoundary(:,1) 4*D0_gasBoundary(:,2)];
D2_gasBoundary=[zeros(lv_gasBoundary,1) zeros(lv_gasBoundary,1) 4*D0_gasBoundary(:,1)];
for j=3:Ng
    D1_gasBoundary=[D1_gasBoundary 2*j*D0_gasBoundary(:,j)+j*D1_gasBoundary(:,j-1)/(j-2)];
    D2_gasBoundary=[D2_gasBoundary 2*j*D1_gasBoundary(:,j)+j*D2_gasBoundary(:,j-1)/(j-2)];
end
    
%% Liquid Interface point    
% Chebyshev matrices are defined for 
% liquid interface (r=R or \tilde{r}_L=+1}

vec_liqInterface = 0;   
rL_interface = cos(vec_liqInterface);   % rL=+1.
    
% The Chebyshev matrices are defined below
D0_liqInterface=[];
for j=0:1:Nl
    D0_liqInterface=[D0_liqInterface cos(j*vec_liqInterface)];
end
lv_liqInterface=length(vec_liqInterface);
D1_liqInterface=[zeros(lv_liqInterface,1) D0_liqInterface(:,1) 4*D0_liqInterface(:,2)];
D2_liqInterface=[zeros(lv_liqInterface,1) zeros(lv_liqInterface,1) 4*D0_liqInterface(:,1)];
for j=3:Nl
    D1_liqInterface=[D1_liqInterface 2*j*D0_liqInterface(:,j)+j*D1_liqInterface(:,j-1)/(j-2)];
    D2_liqInterface=[D2_liqInterface 2*j*D1_liqInterface(:,j)+j*D2_liqInterface(:,j-1)/(j-2)];
end

%% Gas Interface point    
% Chebyshev matrices are defined for 
% gas interface (r=R or \tilde{r}_G=-1}
    
vec_gasInterface = 0;
rG_gasInterface = cos(vec_gasInterface);   % rG=1.
 
% The Chebyshev matrices are defined below 
D0_gasInterface=[];
for j=0:1:Ng
    D0_gasInterface=[D0_gasInterface cos(j*vec_gasInterface)];
end
lv_gasInterface=length(vec_gasInterface);
D1_gasInterface=[zeros(lv_gasInterface,1) D0_gasInterface(:,1) 4*D0_gasInterface(:,2)];
D2_gasInterface=[zeros(lv_gasInterface,1) zeros(lv_gasInterface,1) 4*D0_gasInterface(:,1)];
for j=3:Ng
    D1_gasInterface=[D1_gasInterface 2*j*D0_gasInterface(:,j)+j*D1_gasInterface(:,j-1)/(j-2)];
    D2_gasInterface=[D2_gasInterface 2*j*D1_gasInterface(:,j)+j*D2_gasInterface(:,j-1)/(j-2)];
end

%% %%%%%%%%%%%%%%%%%%%%%%% %%
%%%%% Base velocities %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Liquid

% Base velocity (The equation is non-dimensionalized)
Factor1=kkk*Vr/(Ugstar-qqq);
y1 = (rL+1)/2*R;

y1a= (1+1)/2*R;

Ul = (qqq-Ugstar)*tanh(Factor1*(y1-R))+Ugstar;
Ul = Ul/Ujet;

Ul_int = (qqq-Ugstar)*tanh(Factor1*(y1a-R))+Ugstar;
Ul_int = Ul_int/Ujet;     % interface value

% First derivative of velocity
Ulp = (qqq-Ugstar) * Factor1 * (1./cosh(Factor1 * (y1 - R))).^2;
Ulp = Ulp/Ujet;

Ulp_int = (qqq-Ugstar) * Factor1 * (1./cosh(Factor1 * (y1a - R))).^2;
Ulp_int = Ulp_int/Ujet;   % interface value

% Second derivative of velocity
Ulpp = -2 * (qqq-Ugstar) * Factor1^2 * (1./cosh(Factor1 * (y1 - R))).^2 .* tanh(Factor1 * (y1 - R));
Ulpp = Ulpp/Ujet;

Ulpp_int = -2 * (qqq-Ugstar) * Factor1^2 * (1./cosh(Factor1 * (y1a - R))).^2 .* tanh(Factor1 * (y1a - R));
Ulpp_int = Ulpp_int/Ujet; % interface value

%% Gas

Factor2 = kkk/(Uww-Ugstar); 
y2 = (rG*(1-lr)+(1+lr))/2*R;

y2a= (1*(1-lr)+(1+lr))/2*R;
% Base velocity (The equation is non-dimensionalized)

Ug = (Ugstar-Uww)*tanh(Factor2*(y2-R))+Ugstar;
Ug = Ug/Ujet;

Ug_int = (Ugstar-Uww)*tanh(Factor2*(y2a-R))+Ugstar;
Ug_int = Ug_int/Ujet;     % interface value

% First derivative of velocity
Ugp = (Ugstar-Uww) * Factor2 .* (1./cosh(Factor2 * (y2 - R))).^2;
Ugp = Ugp/Ujet;

Ugp_int = (Ugstar-Uww) * Factor2 .* (1./cosh(Factor2 * (y2a - R))).^2;
Ugp_int = Ugp_int/Ujet;   % interface value

% Second derivative of velocity
Ugpp = -2 * (Ugstar-Uww) * Factor2^2 .* (1./cosh(Factor2 * (y2 - R))).^2 .* tanh(Factor2 * (y2 - R));
Ugpp = Ugpp/Ujet;

Ugpp_int =  -2 * (Ugstar-Uww) * Factor2^2 .* (1./cosh(Factor2 * (y2a - R))).^2 .* tanh(Factor2 * (y2a - R));
Ugpp_int = Ugpp_int/Ujet; % interface value

Uint = Ul_int;  % There is continuity of base velocity at interface (Uint=Ul_int=Ug_int)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%% The following section is solving the equations %%%%%%%
%%%%% for each frequency (alpha) values. %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% System of linear equations in Chebyshev space

for re=1:1
    Rel = Rel_matrix(re,1);
    Wel = Wel_matrix(re,1);
    disp('Jet properties (JFM [2009]):');
    disp(['Density ratio, rho_G/rho_L = ', num2str(Dr)]);
    disp(['Viscosity ratio, mu_G/mu_L = ', num2str(Vr)]);
    disp(['Reynolds number, Re = ', num2str(Rel)]);
    disp(['Weber number, We = ', num2str(Wel)]);
    disp('');
    alphaArray = (ast:mmm:aend)';
    for num=1:length(alphaArray)
      alpha = alphaArray(num,1);    % This is the non-dimensional frequency
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
    %%%%% Construction of matrices (C0,C1,C2) in the eigenvalue %%%%%%%
    %%%%% problem of the form: k^2*C2.a+k*C1.a+C0.a=0 %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
    %%%%% Governing equations for liquid phase %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Continuity
    
    % C0 matrix
    C0_liq_11 = 2*D1_liq + 1/R*(2./(rL+1)).*D0_liq;  % radial velocity part       v'+1/r*v
    C0_liq_13 = 1i*alpha*D0_liq;                 % axial velocity part        i*α*u
    C0_liq_14 = zeros(Nl-1,Nl+1);                % pressure part
    
    % C1 matrix
    C1_liq_11 = zeros(Nl-1,Nl+1);                             
    C1_liq_13 = zeros(Nl-1,Nl+1);                   
    C1_liq_14 = zeros(Nl-1,Nl+1);              
    %% radial momentum
    
    % C0 matrix
    C0_liq_21 = 1/(R^2)*(4/Rel).*D2_liq + 1/(R^2)*(4./(Rel*(rL+1))).*D1_liq - (1/(R^2)*(4./(Rel*((rL+1).^2))) + 1./Rel*alpha^2 + 1i*alpha*Ul).*D0_liq;  
    % radial velocity part      1/Re*v''+1/Re/r*v'-(1/Re/r^2+1/Re*α^2+i*alpha*U)*v
    C0_liq_23 = zeros(Nl-1,Nl+1);                       % axial velocity part
    C0_liq_24 = -1/(R)*2*D1_liq;                              % pressure part   -p'
    
    % C1 matrix
    C1_liq_21 = 1i*D0_liq;          % radial velocity part    i*v*w
    C1_liq_23 = zeros(Nl-1,Nl+1);   % axial velocity part
    C1_liq_24 = zeros(Nl-1,Nl+1);   % pressure part
    
    %% axial momentum
    
    % C0 matrix
    C0_liq_41 = -Ulp.*D0_liq;     % radial velocity part    -U'*v
    C0_liq_43 = 1/(R^2)*(4/Rel)*D2_liq + 1/(R^2)*(4./(Rel*(rL+1))).*D1_liq - 1/Rel*alpha^2*D0_liq - 1i*alpha*Ul.*D0_liq; 
    % axial velocity part 1/Re*u''+1/Re/r*u'-1/Re*α^2*u - i*α*U*u
    C0_liq_44 = -1i*alpha*D0_liq;              % pressure part -i*α*p
    
    % C1 matrix
    C1_liq_41 = zeros(Nl-1,Nl+1);   % radial velocity part
    C1_liq_43 = 1i*D0_liq;          % axial velocity part   i*u*omega
    C1_liq_44 = zeros(Nl-1,Nl+1);   % pressure part
    
    %% Construction of the C0,C1 and C1 matrices for liquid    
    % zero rows at the end is to fill the liquid centerline and interface
    % conditions. (A description of this is given in Slide #57 in the
    % keynote 'Two_Phase_flows_cylindrical_detailed_April_15_2021.key')
    C0_liq = [C0_liq_21 C0_liq_23 C0_liq_24(:,1:end-1);...
        C0_liq_41 C0_liq_43 C0_liq_44(:,1:end-1);...
        C0_liq_11 C0_liq_13 C0_liq_14(:,1:end-1);...
        zeros(5, 3*Nl+2)];
    
    C1_liq = [C1_liq_21 C1_liq_23 C1_liq_24(:,1:end-1);...
        C1_liq_41 C1_liq_43 C1_liq_44(:,1:end-1);...
        C1_liq_11 C1_liq_13 C1_liq_14(:,1:end-1);...
        zeros(5, 3*Nl+2)];
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
    %%%%% Governing equations for gas phase %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Continuity
    
    % C0 matrix
    C0_gas_11 = (2/(1-lr))*D1_gas + (2./(rG*(1-lr)+lr+1)).*D0_gas;	% radial velocity part       v'+1/r*v
    C0_gas_13 = 1i*alpha*D0_gas;                                    % axial velocity part          i*α*u
    C0_gas_14 = zeros(Ng-1,Ng+1);                                   % pressure part
    
    % C1 matrix
    C1_gas_11 = zeros(Ng-1,Ng+1);	% radial velocity part
    C1_gas_13 = zeros(Ng-1,Ng+1);   % axial velocity part
    C1_gas_14 = zeros(Ng-1,Ng+1);   % pressure part
    
    
    %% radial momentum
    
    % C0 matrix
    C0_gas_21 = 1/(R^2)*(Vr/Dr)*(4/(Rel*(1-lr)^2))*D2_gas + 1/(R^2)*(Vr/Dr)*(4./(Rel*(1-lr)*(rG*(1-lr)+lr+1))).*D1_gas -...
        ( 1/(R^2)*(Vr/Dr)*(4./(Rel*((rG*(1-lr)+lr+1).^2))) + (Vr/Dr)/Rel*alpha^2 + 1i*alpha*Ug).*D0_gas;	
    % radial velocity part  1/Re*v''+1/Re/r*v'-(1/Re/r^2+1/Re*α^2+i*alpha*U)*v
    C0_gas_23 = zeros(Ng-1,Ng+1);        % axial velocity part
    C0_gas_24 = -1/(R)*2/((1-lr))/Dr*D1_gas;      % pressure part       -p'
    
    % C1 matrix
    C1_gas_21 = 1i*D0_gas;                 % radial velocity part    i*v*w
    C1_gas_23 = zeros(Ng-1,Ng+1);          % axial velocity part
    C1_gas_24 = zeros(Ng-1,Ng+1);          % pressure part
    
    
    %% axial momentum
    
    % C0 matrix
    C0_gas_41 = -Ugp.*D0_gas;                  % radial velocity part -U'*v
    C0_gas_43 = 1/(R^2)*(Vr/Dr)*(4/(Rel*(1-lr)^2))*D2_gas + 1/(R^2)*(Vr/Dr)/Rel*(4./((1-lr)*(rG*(1-lr)+lr+1))).*D1_gas -...
       (Vr/Dr)/Rel*alpha^2*D0_gas - 1i*alpha*Ug.*D0_gas; 
    % axial velocity part 1/Re*u''+1/Re/r*u'-1/Re*α^2*u- i*α*U*u
    C0_gas_44 = -(1i*alpha)/Dr*D0_gas;           % pressure part -i*α*p
    
    % C1 matrix
    C1_gas_41 = zeros(Ng-1,Ng+1);   % radial velocity part
    C1_gas_43 = 1i*D0_gas;          % axial velocity part i*u*omega
    C1_gas_44 = zeros(Ng-1,Ng+1);   % pressure part
    
    
    %% Construction of the C0,C1 and C1 matrices for gas
    
    % zero rows at the beginning is to fill the gas boundary and interface
    % conditions. 
    
    C0_gas = [zeros(5,3*Ng+2);...
        C0_gas_21 C0_gas_23 C0_gas_24(:,1:end-1);...
        C0_gas_41 C0_gas_43 C0_gas_44(:,1:end-1);...
        C0_gas_11 C0_gas_13 C0_gas_14(:,1:end-1)];
    
    C1_gas = [zeros(5,3*Ng+2);...
        C1_gas_21 C1_gas_23 C1_gas_24(:,1:end-1);...
        C1_gas_41 C1_gas_43 C1_gas_44(:,1:end-1);...
        C1_gas_11 C1_gas_13 C1_gas_14(:,1:end-1)];
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
    %%%%% Boundary condition at liquid centerline %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % radial velocity

    tempC0 = [D0_liqCenter zeros(1,Nl+1) zeros(1,Nl)];    %  v=0
    tempC1 = [zeros(1,Nl+1) zeros(1,Nl+1) zeros(1,Nl)];
    
    C0_liq(3*Nl-2,:) = tempC0;
    C1_liq(3*Nl-2,:) = tempC1;
    
    % axial velocity
    
    tempC0 = [zeros(1,Nl+1) D1_liqCenter zeros(1,Nl)];   % u'=0
    tempC1 = [zeros(1,Nl+1) zeros(1,Nl+1) zeros(1,Nl)];
  
    C0_liq(3*Nl-1,:) = tempC0;
    C1_liq(3*Nl-1,:) = tempC1;
    
    % pressure

    tempC0 = [zeros(1,Nl+1) zeros(1,Nl+1) D1_liqCenter(:,1:end-1)]; % p'=0
    tempC1 = [zeros(1,Nl+1) zeros(1,Nl+1) zeros(1,Nl)];
   
    C0_liq(3*Nl,:) = tempC0;
    C1_liq(3*Nl,:) = tempC1;
   
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
    %%%%% Boundary condition at gas boundary %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % radial velocity

    tempC0 = [D0_gasBoundary zeros(1,Ng+1) zeros(1,Ng)];    %  v=0
    tempC1 = [zeros(1,Ng+1)  zeros(1,Ng+1) zeros(1,Ng)];
    
    C0_gas(5,:) = tempC0;
    C1_gas(5,:) = tempC1;    
    
    % axial velocity
    
    tempC0 = [zeros(1,Ng+1) D1_gasBoundary zeros(1,Ng)];   %  u'=0
    tempC1 = [zeros(1,Ng+1) zeros(1,Ng+1) zeros(1,Ng)];  
   
    C0_gas(4,:) = tempC0;
    C1_gas(4,:) = tempC1;
    
    % pressure
    
    tempC0 = [zeros(1,Ng+1) zeros(1,Ng+1) D1_gasBoundary(:,1:end-1)];  %  p'=0
    tempC1 = [zeros(1,Ng+1) zeros(1,Ng+1) zeros(1,Ng)];
   
    C0_gas(3,:) = tempC0;
    C1_gas(3,:) = tempC1;
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
    %%%%% Combined matrices for the two phases %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % This is excluding the interface conditions.
    
    C0 = [C0_liq zeros(3*Nl+2,3*Ng+2); zeros(3*Ng+2,3*Nl+2) C0_gas];
    C1 = [C1_liq zeros(3*Nl+2,3*Ng+2); zeros(3*Ng+2,3*Nl+2) C1_gas];
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
    %%%%% Interface conditions %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
    %% Jump in Normal stress condition  压力平衡

    tempC0_liq = [-1/(R)*4/Rel*D1_liqInterface 1/Wel*(alpha^2-1)*(1/(Ulp_int-Ugp_int))*D0_liqInterface D0_liqInterface(:,1:end-1)];
    tempC1_liq = [zeros(1,Nl+1) zeros(1,Nl+1) zeros(1,Nl)];
    % 
    tempC0_gas = [1/(R)*4*Vr/(Rel*(1-lr))*D1_gasInterface 1/Wel*(-alpha^2+1)*(1/(Ulp_int-Ugp_int))*D0_gasInterface -D0_gasInterface(:,1:end-1)];
    tempC1_gas = [zeros(1,Ng+1) zeros(1,Ng+1) zeros(1,Ng)];

     
    C0(3*Nl+1,:) = [tempC0_liq tempC0_gas];
    C1(3*Nl+1,:) = [tempC1_liq tempC1_gas];
    
    %% Continuity of Shear stress condition-1  剪切力平衡  
    


     tempC0_liq = [-(alpha^2*Uint)*D0_liqInterface 1/(R)*1i*alpha*Uint*2*D1_liqInterface zeros(1,Nl)];
     tempC1_liq = [alpha*D0_liqInterface -1/(R)*1i*2*D1_liqInterface zeros(1,Nl)];
     % 
     tempC0_gas = [Vr*alpha^2*Uint*D0_gasInterface -1/(R)*1i*alpha*Vr*Uint*2/(1-lr)*D1_gasInterface zeros(1,Ng)];   %%%%%???
     tempC1_gas = [-alpha*Vr*1i*D0_gasInterface 1/(R)*1i*Vr*2/(1-lr)*D1_gasInterface zeros(1,Ng)];


   
    C0(3*Nl+2,:) = [tempC0_liq tempC0_gas];
    C1(3*Nl+2,:) = [tempC1_liq tempC1_gas];
    %% Continuity of radial velocity condition    径向速度 -v1+v2=0
    
    tempC0_liq = [-D0_liqInterface  zeros(1,Nl+1) zeros(1,Nl)];
    tempC1_liq = [zeros(1,Nl+1)  zeros(1,Nl+1) zeros(1,Nl)];
    
    tempC0_gas = [D0_gasInterface zeros(1,Ng+1) zeros(1,Ng)];
    tempC1_gas = [zeros(1,Ng+1) zeros(1,Ng+1) zeros(1,Ng)];

    C0(3*Nl+3,:) = [tempC0_liq tempC0_gas];
    C1(3*Nl+3,:) = [tempC1_liq tempC1_gas];
    %% Continuity of axial velocity condition 轴向速度平衡 u1-u2+（U1'-U2')*nt=0
    tempC0_liq = [(Ulp_int-Ugp_int)*D0_liqInterface alpha*1i*Uint*D0_liqInterface zeros(1,Nl)];
    tempC1_liq = [zeros(1,Nl+1) -1i*D0_liqInterface zeros(1,Nl)];
    
    tempC0_gas = [zeros(1,Ng+1) -alpha*1i*Uint*D0_gasInterface zeros(1,Ng)];
    tempC1_gas = [zeros(1,Ng+1) 1i*D0_gasInterface zeros(1,Ng)];

    C0(3*Nl+4,:) = [tempC0_liq tempC0_gas];
    C1(3*Nl+4,:) = [tempC1_liq tempC1_gas];
    
 
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
    %%%%% Solving for the eigenvalues %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    % The non-linear eigenvalue system k^2*C2.a+k*C1.a+C0.a=0 is linearized
    % iw*C1*F+C0*F=0;  
    % C0*F=-w*C1*F;
    
    % 改进后的特征值追踪逻辑

    % 1. 动态选择 eigs 的搜索起点（Target Value）
    % 初始时设为 0（原代码的设置），后续使用上一步追踪到的特征值。
    if num == 1
        % 第一次计算：使用 0 作为搜索起点
         alpha_target = ast;
    else
         % 后续计算：使用上一次追踪到的特征值作为搜索起点，提高追踪精度和速度
         % 注意：kFinal(num-1, 1) 是上一步存储的特征值 omega
        alpha_target = kFinal(num-1, 1);
     end
     C1 = C1 + eye(size(C1)) * 1e-9;
     % 2. 求解广义特征值问题
     % 求解 50 个最接近 alpha_target 的特征值
     [F_curr, k_curr_diag] = eigs(C0, -C1, 50, alpha_target); 

     
     k_curr = diag(k_curr_diag); % 提取特征值 omega
     [~, M_eigs] = size(F_curr); % M_eigs 是求解的特征值数量

    
         % 选取增长率（虚部）最大的模态
         [M, I] = max(imag(k_curr)); 
    
         F_tracked = F_curr(:, I);
         k_tracked_val = k_curr(I);

      

       % 9. 存储最终结果
       kFinal(num, 1) = k_tracked_val;      % 存储复特征值 omega
       Ffinal(:,num) = F_tracked;           % 存储追踪到的特征向量
       
       
       GrowthRate(num,1) = imag(k_tracked_val); % 增长率 (虚部)
     
end

%% %%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%% Dispersion plots %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure();
plot(alphaArray, GrowthRate,'o--')

xlabel('$Wavenumber$','Interpreter','Latex')
ylabel('$Growthrate$','Interpreter','Latex')
title({
    ['JFM2009']
    ['Re=' num2str(Rel) ' ,We=' num2str(Wel) ' ,K=',num2str(kkk) ' ,Us=',num2str(Ugstar)]
    },'FontSize',20);
grid on
box on

    end
