%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Layers Structure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first column is the material used from the "library"
% second column is the length of the layer in nm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%substrate=GaAs;      % Important for the Strain model (GaAs, InP, InAs, GaSb)
%M=[
%GaAs      10
%InGaAs20  10
%GaAs      10
%];

%substrate=GaAs;      % Important for the Strain model (GaAs, InP, InAs, GaSb)
%M=[
%AlGaAs40      5
%GaAs          5
%AlGaAs40      5
%
%];

% Sirtori, PRB, 50, 8663 (1994)
%substrate=InP;      % Important for the Strain model (GaAs, InP, InAs, GaSb)
%M=[
%AlInAs   6
%InGaAs   5.2
%AlInAs   6
%];

% Sirtori, PRB, 50, 8663 (1994)
%substrate=InP;      % Important for the Strain model (GaAs, InP, InAs, GaSb)
%M=[
%AlInAs   6
%InGaAs   5.9
%AlInAs   1.3
%InGaAs   2.4
%AlInAs   6
%];

% Sirtori, PRB, 50, 8663 (1994)
substrate=InP;      % Important for the Strain model (GaAs, InP, InAs, GaSb)
M=[
AlInAs   6
InGaAs   4.6
AlInAs   1.0
InGaAs   2.0
AlInAs   1.0
InGaAs   1.9
AlInAs   6
];

%substrate=InP;      % Important for the Strain model (GaAs, InP, InAs, GaSb)
%M=[
%GaAsSb   10
%InGaAs   10
%GaAsSb   10
%];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
