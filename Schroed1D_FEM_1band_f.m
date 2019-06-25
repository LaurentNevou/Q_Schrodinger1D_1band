%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Schrodinger solver on uniform grid with m(z)!!! %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[E,psi]=Schroed1D_FEM_1band_f(z,V0,Mass,n)

% One problem here is the operator (d/z)*(m(z))*(d/z) inside the Shrodinger equation
% Without any modification, the results are not accurate and show discontinuity 
% at the material interface due to the difference of masses
% One very nice approach is to use the mid-point mass:

% Paul Harrisson
% Quantum Wells, Wires and Dots.
% 4th edition (2016),
% chap 3 : "Numerical Solutions"
% 3.13: "Extention to variable effective mass"
% page 102

% Xunpeng Ma, Kangwen Li, Zuyin Zhang, Haifeng Hu, Qing Wang, Xin Wei,and Guofeng Song
% "Two-band finite difference method for the bandstructure calculation with nonparabolicity effects in quantum cascade lasers"
% JAP, 114, 063101 (2013)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=6.62606896E-34;               %% Planck constant J.s
hbar=h/(2*pi);
e=1.602176487E-19;              %% charge de l electron Coulomb
m0=9.10938188E-31;              %% electron mass kg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nz=length(z);
dz = z(2)-z(1);

shift=min(V0);
V0=V0-shift;
V0(1) = 2;
V0(end) = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Building of the operators %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Massp = [         (Mass(1:end-1) + Mass(2:end)) / 2   Mass(end) ];
Massm = [ Mass(1) (Mass(1:end-1) + Mass(2:end)) / 2              ];

b = (1./Massp + 1./Massm) .* ones(1,Nz) ;
a = 1./Massm(2:end) .* ones(1,Nz-1) ;
c = 1./Massm(2:end) .* ones(1,Nz-1) ;

DZ2 = (-1)*diag(b)  +  (1)*diag(a,-1)  +  (1)*diag(c,+1) ;

DZ2 = DZ2 / dz^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Building of the Hamiltonien %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H = (-hbar^2/(2*m0)) * DZ2 + diag(V0*e);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Diagonalisation of the Hamiltonien %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H=sparse(H);
[psi,Energy] = eigs(H,n,'SM');
E = diag(Energy)/e + shift;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Normalization of the Wavefunction %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:n
  psi(:,i)=psi(:,i)/sqrt(trapz(z',abs(psi(:,i)).^2));  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% here is a small patch due to differences between Octave and Matlab
% Matlab order the eigen values while Octave reverse it

if E(1)>E(2)
  psi=psi(:,end:-1:1);
  E=E(end:-1:1);
end

end