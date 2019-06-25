%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% last update 25June2019, lne %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program solves the Schrodinger equation with m(z) using 2 different algorithms:
% FDM (Finite Difference Method) or scaning/shooting methods
%
% -> A strain model is included. It basically shifts the conduction band edge
% It does not influence the mass since it is a 1 band model, m=m(z). 
% The strain is mainly interesting for InGaAs/GaAs heterostructures

% -> Si/Ge, II-VI and cubic nitride material parameters are available but should
% be grabt in the "Library.m" file
% -> Wurtzite parameter are also availables but the code isn't optimized for it.
% In the "Library.m" file, the WZ table must be open and the meaningfull parameters
% must be taken. Also, the electric field has to be handle...
% -> Additionnal material can be added in the "materialDB_ZB.csv" file

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=6.62606896E-34;               %% Planck constant J.s
hbar=h/(2*pi);
e=1.602176487E-19;              %% charge de l electron Coulomb
m0=9.10938188E-31;              %% electron mass kg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0 for turn off
% 1 for turn on

FE_1band_Method    = 1;       % Diagonalization of the Hamiltonian (FEM)
Shooting_1band_Method = 0;    % Scanning in Energy method (Euler)

StrainModel  = 0;             % Activate Strain model (Shift the CBand only)
PrintResults = 1;             % Switch to print or not the dipoles on the shell

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=3;                      % number of solution asked
ScF=0.1;                  % scaling factor to plot the wave function [Without Dimension]
dz=0.5E-10;               % resolution of the grid [m]
F0=0;%-60e6;%-9.          % Electric field [Volt/meter]
T=300;                    % Temperature [Kelvin], react on the band gap Eg only

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Library;                  % load material parameter DB from "materialDB_ZB.csv"
ExtractParameters;        % extract parameter from the Library
TernaryAlloy;             % compute the ternary alloy    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% import the layer structure file %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first column is the material used from the "library"
% second column is the length of the layer in nm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

input_file;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% NOTHING TO CHANGE ANYMORE !!! %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Grabbing the parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zt  = M(:,end)*1E-9;      % conversion of the length from Angstrom to meter

Egt = M(:,idx_Eg6c) - (M(:,idx_alphaG)*T^2) ./ (T+M(:,idx_betaG));   %Eg = Eg0 - (a*T.^2)./(T + b);
VBOt= M(:,idx_VBO);
met = M(:,idx_me);

if StrainModel == 1

at  = M(:,idx_a);           % lattice parameter
act = M(:,idx_ac);          % Conduction band strain offset parameter
c11t = M(:,idx_c11);        % Conduction band strain offset parameter
c12t = M(:,idx_c12);        % Conduction band strain offset parameter

a0 = substrate(idx_a);
exx =  (a0-at)/a0;
eyy =  exx;
ezz = -2*c12t./c11t.*exx;

DCBOt = -abs(act).*(exx+eyy+ezz);

CBOt= Egt + VBOt + DCBOt ;  % CBO form band gap difference and temperature

else

CBOt= Egt + VBOt ;  % CBO form band gap difference and temperature

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Discretisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% here, I descretize the grid z, the potential V0 and the mass me

z(1)=0;  V0(1)=CBOt(1); Eg(1)=Egt(1); me(1)=met(1);

for i=1:length(zt)
    zv  =  (z(end)+dz) : dz : (z(end) + zt(i)) ;
    z   = [ z  zv ];
    V0  = [ V0     ones(size(zv)) * CBOt(i)  ];
    me  = [ me     ones(size(zv)) * met(i)   ];
end

V0=V0-min(V0);
V0=(F0*z)+V0;              % adding the electric field to the potential

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Selection of the model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if FE_1band_Method==0 && Shooting_1band_Method==0
  display('ERROR: Choose one of the model, FEM or Shooting method')
  return
end

s{1}='\fontsize{20}\color{blue}Potential';
s{2}='\fontsize{20}\color{blue}----------';
j=0;FS=20;
c=[
1 0 0
0 1 0
0 0 1
];

if FE_1band_Method==1
    tic
    j=j+1;
    Mass=me; % m=m(z)
    [E{j},psi{j}] = Schroed1D_FEM_1band_f(z,V0,Mass,n);  % m = m(z)
    s{j+2}=strcat('\fontsize{',num2str(FS),'}\color[rgb]{',num2str(c(j,:)),'}-> FDM-1band: m=m(z)');
    display(strcat('-> FDM 1band method =',num2str(toc),'sec'))
end

if Shooting_1band_Method==1
    dE=0.02;
    precision=1e-5;
    tic
    j=j+1;
    Mass=me; % m=m(z)
    s{j+2}=strcat('\fontsize{',num2str(FS),'}\color[rgb]{',num2str(c(j,:)),'}-> Shooting-1band: m=m(z)');
    [E{j},psi{j}] = Schroed1D_Euler_1band_f(z,V0,Mass,n,dE,precision);
    display(strcat('-> Shooting 1band method =',num2str(toc),'sec'))      
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Display Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=1:length(E)
  for i=1:length(E{k})
      PSI{k}(:,i)=abs(psi{k}(:,i)).^2/max(abs(psi{k}(:,i)).^2)*ScF + E{k}(i); % normalisation for the plotting
  end
end

for k=1:length(E)
  for i=1:length(E{k})
    for j=1:length(E{k})
      if j>i
        z_dipole{k}(i,j) = abs(  trapz( z , psi{k}(:,i).*z'.*psi{k}(:,j) )  );
        f_dipole{k}(i,j) = 2*m0/hbar^2 * ( E{k}(j)-E{k}(i) )* e * z_dipole{k}(i,j)^2 ;
        % Take care! Some people use meff inside the oscillator strenght f
        % Actually, meff has sens in an infinite QW because there is a single mass value
        % but not in multi-QW structure with various materials
        % https://www.nextnano.com/nextnano3/tutorial/1Dtutorial_IntrabandTransitions.htm
        % https://www.nextnano.com/nextnano3/tutorial/1Dtutorial_InGaAs_MQWs.htm
      end
    end
  end
end

if PrintResults == 1

display('============================================')
display('Results:')
display('============================================')

for k=1:length(E)

if k>1
  display('============================================')
end

for i=1:length(E{k})
  for j=1:length(E{k})
    if j>i
      
        display(strcat(...
        'E',num2str(i),'-',num2str(j),' = ',num2str( E{k}(j)-E{k}(i),'%.3f' ),'eV;   z'...
        ,num2str(i),'-',num2str(j),' = ',num2str( z_dipole{k}(i,j)*1e9,'%.3f' ),'nm;   f'...
        ,num2str(i),'-',num2str(j),' = ',num2str( f_dipole{k}(i,j),'%.3f' ) ...
        )  )
        
    end
  end
end
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%figure('position',[-1500 100 1000 700]);
figure('position',[10 100 1000 700]);
subplot(1,1,1,'fontsize',15)
hold on;%grid on;

yscale=[min(V0)-0.1 max(V0)+0.1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(z*1e9,V0, 'b-','linewidth',2)

for j=1:length(E)
  for i=1:length(E{j})
    
      plot(z*1e9,PSI{j}(:,i),'color',c(j,:),'linewidth',2)
      
  end
end

xlabel('z (nm)');
ylabel('Energy (eV)');

ylim(yscale)

text(0.5,(yscale(2)-yscale(1))*0.1+yscale(1),s);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
