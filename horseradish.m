% horseradish:   Calculate EDNMR spectrum for arbitrary spin system
%
%   [nu,spec] = horseradish(Sys,Exp)
%   ... = horseradish(Sys,Exp,Opt)
%
%   Input:
%     Sys      spin system structure
%     Exp      experiment structure
%       Exp.Field        magnetic field, in mT
%       Exp.Range        mw offset range [min max], in MHz
%       Exp.nPoints      number of points
%       Exp.ExciteWidth  excitation bandwidth, in MHz
%       Exp.Q            Q0 of the cavity
%       Exp.tHTA         HTA pulse length, in us
%       Exp.nu1          nu1, in MHz, normalized for spin1/2 with g=gfree
%       Exp.Tm           decay time of electron coherences, in us
%       Exp.Temperature
%     Opt      options structure
%       Opt.Symmetry     symmetry of spin system
%       Opt.nKnots       number of knots for the orientation grid
%       Opt.Threshold.Probe    cutoff for transition selection
%       Opt.Threshold.Pump     cutoff if HTA has an effect on

function [nu,spec] = horseradish(varargin)

if (nargin==0), help(mfilename); return; end
[Sys,Exp,Opt]=ParseInputs(varargin{:});

% Create Orientational Grid
%%% hey nino %%% - you should add possibility of single crystals
[phi,theta,PowderWeight] = sphgrid(Opt.Symmetry,Opt.nKnots);
nOrientations = numel(PowderWeight);

% Pre-calculate spin Hamiltonian components in molecular frame
[Hzf,GxM,GyM,GzM] = sham(Sys);
zeemanop = @(v) v(1)*GxM + v(2)*GyM + v(3)*GzM;

% initialize output
%%% hey nino %%% - think about possible outputs as per orientation
[nu,spec]=makespec(Exp.Range,Exp.nPoints,Exp.Range,[0 0]);

% loop over orientations
for iOrient=1:nOrientations
    
    Exp.CrystalOrientation = [phi(iOrient) theta(iOrient) 0];
    % calculate energy levels and probabilites between them
    [EnergyLevels, Probabilities]=SolveEigenProblem(Hzf,zeemanop,Exp,Opt);
    
    % get populations
    %%% hey nino %%% - make user-supplied starting states possible
    Populations=diag(  sigeq(diag(EnergyLevels),Exp.Temperature,'pol')  );
    
    % transition frequencies between levels in MHz
    TransFreqs = abs(EnergyLevels-EnergyLevels');
    TransOffsets = (TransFreqs-Exp.mwFreq);
    
    %calculate probe weight of each transition -  ensures orientation
    %selection
    ProbeWeights = Probabilities.*exp(-2*((TransFreqs-Exp.mwFreq)/Exp.ExciteWidth).^2);
    ProbeWeights(ProbeWeights<Opt.Threshold.Probe)=0; %apply cutoff
    
    if ~any(ProbeWeights)
        continue;
    end
    
    % Calculate W1(HTA) at EDNMR position using experimental loaded Q-value
    w1 = 2*pi*Exp.nu1*sqrt(1./(1+(TransOffsets*4*Exp.Q/(Exp.mwFreq)).^2/4));
    PumpWeights=(1-Bloch(w1,Exp.tHTA,Probabilities,Exp.Tm))/2; %between 0, and 1 (inversion)
    PumpWeights(PumpWeights<Opt.Threshold.Pump)=0;  % apply cutoff
    PumpWeights(TransOffsets<Exp.Range(1) | TransOffsets>Exp.Range(2))=0; % only hit when within range
    
    [Amp, Pos]=EDNMR_PeakList(TransFreqs,ProbeWeights,PumpWeights,Populations);
    
    %add up the spectrum from this powder orientation
    if ~isempty(Pos)
        spec=spec+PowderWeight(iOrient)*makespec(Exp.Range,Exp.nPoints,Pos,Amp);
    end
    
end %end of powder loop

dnu = nu(2)-nu(1);
spec = convspec(spec,dnu,Sys.lwEndor);

end % end of main function

function [Sys,Exp,Opt]=ParseInputs(varargin)

% Guard against wrong number of input or output arguments.
if (nargin<1), error('Please supply a spin system as first parameter.'); end
if (nargin<2), error('Please supply experimental parameters as second input argument.'); end
if (nargin>3), error('Too many input arguments, the maximum is three.'); end

if (nargout>4), error('Too many output arguments.'); end

% Put vargin into structures, check, and add default values
Sys=varargin{1};
Exp=varargin{2};
%%% hey nino %%% - you should add some more Exp-checks here
if ~isfield(Exp,'Temperature') Exp.Temperature=300; end
Exp.mwFreq=Exp.mwFreq*1e3; %GHz to MHz

if nargin==3
    Opt=varargin{3};
else
    Opt=struct;
end
if ~isfield(Opt,'Symmetry') Opt.Symmetry=symm(Sys); end
if ~isfield(Opt,'nKnots') Opt.nKnots=31; end
if ~isfield(Opt,'Symmetry') Opt.Symmetry=symm(Sys); end
if ~isfield(Opt,'Threshold')
    Opt.Threshold=struct;
end
if ~isfield(Opt,'Threshold.Probe')  Opt.Threshold.Probe=1e-4; end
if ~isfield(Opt,'Threshold.Pump')  Opt.Threshold.Pump=1e-4; end


end

function [EnergyLevels, Probabilities]=SolveEigenProblem(Hzf,zeemanop,Exp,Opt)

% Set up lab axes for this orientation
[xL,yL,zL] = erot(Exp.CrystalOrientation,'rows');
% zL = z laboratoy axis: external static field
% xL = x laboratory axis: B1 excitation field
% yL = y laboratory vector: needed for intensity integration over third Euler angle

% Set up spin Hamiltonian and calculate energy levels and transition
% probabilities (normalized to 1 for a spin-1/2 with g=gfree)
H0 = Hzf + Exp.Field*zeemanop(zL);
[U,Hd] = eig(H0);
EnergyLevels = diag(Hd); %Energy levels

gamma=gfree*bmagn/planck*1e-9; %MHz/mT
Probabilities = 2/gamma^2*(abs(U'*zeemanop(xL)*U).^2+abs(U'*zeemanop(yL)*U).^2);
%%% hey nino %%%- is this the right way to integrate over the third angle?

end

function [Mz] = Bloch(w1,t,P,Tm)

%for now assuming  [0 0 1] as startingvector

fac = exp(-t/(2*Tm));
Ch = cosh(t*sqrt(1-4*P.*Tm.^2.*w1.^2)/(2*Tm));
Sh = sinh(t*sqrt(1-4*P.*Tm.^2.*w1.^2)/(2*Tm));
Mz = fac.*(Ch + Sh./(sqrt(1-4*P.*Tm.^2.*w1.^2)));

end

function [Amp, Pos]=EDNMR_PeakList(TransFreqs,ProbeWeights,PumpWeights,Populations)

%initalize output
Pos=[];
Amp=[];

%get flip anlge from weight
CosPump=-(2*PumpWeights-1);

%abbreviations
Pops=Populations;

%loop over all probe transition
for i_probe=1:size(TransFreqs,1)
    for j_probe=(i_probe+1):size(TransFreqs,2)
        
        if ProbeWeights(i_probe,j_probe)
            
            %loop over all pump transition
            for i_pump = 1:size(TransFreqs,1)
                for j_pump=(i_pump+1):size(TransFreqs,2)
                    
                    %only transitions connected to the probe and affected
                    %by HTA are considered
                    if numel(unique([i_probe j_probe i_pump j_pump]))==3 && PumpWeights(i_pump,j_pump)
                        
                        Pos=[Pos TransFreqs(i_pump,j_pump)-TransFreqs(i_probe,j_probe)]; %postion of pumped transition wrt to observed one(not! the resonator. is this correct?)
                        
                        %calculate polarzation before and after pump
                        %pulse
                        PumpPops = Pops;
                        PumpPops(i_pump)=CosPump(i_pump,j_pump)*Pops(i_pump)+(1-CosPump(i_pump,j_pump))*Pops(j_pump);
                        PumpPops(j_pump)=CosPump(i_pump,j_pump)*Pops(j_pump)+(1-CosPump(i_pump,j_pump))*Pops(i_pump);
                        
                        PolDif = (PumpPops(j_probe)-PumpPops(i_probe))-(Pops(j_probe)-Pops(i_probe));
                        
                        %calculate peak amplitude
                        Amp=[Amp PolDif*ProbeWeights(i_probe,j_probe)];
                    end
                end
            end
        end
    end
end
end