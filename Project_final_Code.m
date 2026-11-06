clc;
clear ;
close all;
%%
disp("    ==============================================================================")
disp("                                     Team MATLITE  ")
disp("    ==============================================================================")
%% Define units
global m nm kg sec J eV mass_e h_bar well_width check
m = 1;
kg = 1;
sec = 1;
J = 1;
nm = 1E-9 *(m);
%%

eV = 1.60218E-19 *(J);

mass_e = 9.10938E-31 *(kg);
h_bar = (6.62607E-34)/(2*pi) *(J*sec);



%% Taking inputs

fprintf('    Potential Well Type: ')
fprintf( '\n\n        i.Input 1 for Finite \n       ii.Input 2 for Ininite\n\n');
fprintf('       >>');
num = input('');

%%
well_width=0;
high_p=0;
x_rag=0;
dmin=0;

check = 0;

while(is_inside(1*nm,20*nm,well_width))
if(check) 
    disp('Invalid input')
end
fprintf('\r    Potential Well Width(between 1-20 (nm)): ');
well_width = input('')*(nm);
check = is_inside(1*nm,20*nm,well_width);
end

check = 0;

while(is_inside(0.01*eV,1*eV,high_p))
if(check)
    disp('Invalid input')
end
fprintf('\r    Higher potential (between 0.01-1 (eV)): ');
high_p = input('')*(eV);
check = is_inside(0.01*eV,1*eV,high_p);
end

check = 0;

while(is_inside(10*nm,50*nm,x_rag))
if(check)
    disp('Invalid input')
end
fprintf('\r    Full potential well length  (between 10-50 (nm)): ')
x_rag = input('')*(nm);
check = is_inside(10*nm,50*nm,x_rag);
end

check = 0;

while(is_inside(0.001*nm,0.1*nm,dmin))
if(check) 
    disp('Invalid input')
end
fprintf('\r    Lowest difference (between 0.001-0.1 (nm)): ')
dmin = input('')*(nm);
check = is_inside(0.001*nm,0.1*nm,dmin);
end

%% For Finite Well
if num ==1
   
%% Determinate simulation parameters
dx = dmin;
xRange = x_rag;
qWell_wide = well_width;
V_upper = high_p;
V_lower = 0.00 *(eV);

%% Preallocating variables and functions
x = (-xRange/2:dx:xRange/2); 

%% constructure the quantum well
V = zeros(size(x));

V(x<-qWell_wide/2) = V_upper;
V(x>qWell_wide/2) = V_upper;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V(x>=0 & x<=qWell_wide/2) = V_lower;
V(x<=0 & x>=-qWell_wide/2) = V_lower;

figure(1), plot(x/nm, V/eV,'LineWidth',5)
xlim([min(x/nm),max(x/nm)]), ylim([min(V/eV),max(1.4*V/eV)]),grid on
 title('FiniteWell')
xlabel('x [nm]'),ylabel('V (eV)')

%% construct the Hamiltonian operator
H0 = diag((h_bar^2/(mass_e*dx^2)) + V);
H1 = diag(ones(1,length(x)-1)*(-h_bar^2/(2*mass_e*dx^2)),1);
H2 = diag(ones(1,length(x)-1)*(-h_bar^2/(2*mass_e*dx^2)),-1);
H = H0+H1+H2;

%% Calculate field propagation
[psi,En] = eig(H);

%%
plotModes = sum(diag(En)>V_lower & diag(En)<V_upper);
figure(1)
for plotMode = 1:plotModes
    rescale = ((En(1,1)-V_lower))/(max(psi(:,plotMode))-min(psi(:,plotMode)));
    psi_rescal_plot = psi(:,plotMode)*rescale + En(plotMode,plotMode);
    hold on, plot(x/nm, psi_rescal_plot/eV,'LineWidth',3)
    hold on, plot(x/nm, En(plotMode, plotMode)*ones(size(x))/eV,'LineWidth',1.2)
    
end


end


%% For Infinite Well

if num ==2
 %% Determinate simulation parameters
dx = dmin;
xRange = x_rag;
qWell_wide = well_width;
V_upper = high_p;
V_lower = 0.00 *(eV);

%% Preallocating variables and functions

x = (-xRange/4:dx:xRange/4);
y = (-xRange:dx:xRange);

%% constructure the quantum well
V = zeros(size(x));

V(y<-qWell_wide/2) = V_upper;
V(y>qWell_wide/2) = V_upper;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V(y>0 & y<qWell_wide) = V_lower;
V(y<0 & y>-qWell_wide) = V_lower;

figure(1), plot(y/nm, V/eV,'LineWidth',4)
xlim([min(y/nm),max(y/nm)]), ylim([min(V/eV),max(1.1*V/eV)]),grid on
title('InfiniteWell')
xlabel('x [nm]'),ylabel('V (eV)')

%% construct the Hamiltonian operator
H0 = eye(length(x),length(x))*(-2)*(-h_bar^2/(2*mass_e*dx^2));
H1 = diag(ones(1,length(x)-1)*(-h_bar^2/(2*mass_e*dx^2)),1);
H2 = diag(ones(1,length(x)-1)*(-h_bar^2/(2*mass_e*dx^2)),-1);
H = H0+H1+H2;

%% Calculate field propagation
[psi,En] = eig(H);

%%
plotModes = sum(diag(En)>V_lower & diag(En)<V_upper);
figure(1)
for plotMode = 1:plotModes
    rescale = ((En(1,1)-V_lower))/(max(psi(:,plotMode))-min(psi(:,plotMode)));
    psi_rescal_plot = psi(:,plotMode)*rescale + En(plotMode,plotMode);
    E = En(plotMode, plotMode)*ones(size(x));
    hold on, plot(x/nm, psi_rescal_plot/eV,'LineWidth',3)
    hold on, plot(x/nm, E/eV,'LineWidth',1.2)        
end
end

%%
disp('')
disp('    Team Members-')
disp('')
disp('        1.Ilias Ahamd (2006068)')
disp('        2.Nahian Tasnim (2006098)');

disp('')
disp('    Thank you')

%% checking Function
function y = is_inside(a,b,c)
  if( a<=c && b>=c ) 
      y = 0;
  
  else
      y = 1;
end
end
