% cheesy:   CHEESY-detected nmr spectrum for arbitrary spin system
%
%   [nu,spec] = cheesy(Sys,Exp)
%   ... = cheesy(Sys,Exp,Opt)
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
%       Exp.Temperature  Temperature in K for inital density matrix
%     Opt      options structure
%       Opt.Symmetry     symmetry of spin system
%       Opt.nKnots       number of knots for the orientation grid
%       Opt.Output       for Powders: 'summed': sum of isotopologues, else
%                        'separate'
%       Opt.Threshold.Probe    cutoff for transition selection
%       Opt.Threshold.Pump     cutoff if HTA has an effect on
%       Opt.Threshold.Iso      cutoff for inclusion of isotopologues

function [nu,spec] = cheesy(varargin)

if (nargin==0), help(mfilename); return; end
[Sys,Exp,Opt]=ParseInputs(varargin{:});

summedOutput=0;
if strcmp(Opt.Output,'summed')
    summedOutput=1;
end

if ~isfield(Sys,'singleiso') || ~Sys.singleiso
    
    
    SysList = isotopologues(Sys,Opt.Threshold.Iso);
    nIsotopologues = numel(SysList);
    
    PowderSimulation = ~isfield(Exp,'CrystalOrientation')  || isempty(Exp.CrystalOrientation);
    
    if ~PowderSimulation && ~summedOutput && nIsotopologues<2
        warning('There is no support for separate output of several crystal orientations without isotopologues. Ignored.')
    end
    
    appendSpectra = (PowderSimulation && ~summedOutput) || (~PowderSimulation && ~summedOutput && nIsotopologues>1);
    if appendSpectra
        spec = [];
    else
        spec = 0;
    end
    
    % Loop over all components and isotopologues
    for iIsotopologue = 1:nIsotopologues
        
        % Simulate single-isotopologue spectrum
        Sys_ = SysList(iIsotopologue);
        Sys_.singleiso = true;
        [nu,spec_] = cheesy(Sys_,Exp,Opt);
        
        % Accumulate or append spectra
        if appendSpectra
            spec = [spec; spec_*Sys_.weight];
        else
            spec = spec + spec_*Sys_.weight;
        end
        
    end
    return
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Simulation of single isotopologue
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create Orientational Grid
%%% hey nino %%% - this is not solved nicely, since Exp.CO writes into phi
%%% and theta, then back further down, doesnt make any sense. but now
%%% needed because i use several functions
if ~isfield(Exp,'CrystalOrientation') || isempty(Exp.CrystalOrientation)
    Opt.Symmetry=symm(Sys);
    [phi,theta,PowderWeight] = sphgrid(Opt.Symmetry,Opt.nKnots);
else
    phi=Exp.CrystalOrientation(:,1);
    theta=Exp.CrystalOrientation(:,2);
    %several crystal orientations are weighted the same for now
    PowderWeight=ones(size(Exp.CrystalOrientation,1),1); 
end
nOrientations = numel(PowderWeight);

% Pre-calculate spin Hamiltonian components in molecular frame
[Hzf,GxM,GyM,GzM] = sham(Sys);
zeemanop = @(v) v(1)*GxM + v(2)*GyM + v(3)*GzM;

% initialize output
%%% hey nino %%% - think about possible outputs as per orientation
[nu,spec]=makespec(Exp.Range,Exp.nPoints,Exp.Range+[eps -eps],[0 0]);

% loop over orientations
parfor iOrient=1:nOrientations
    
    Expl=Exp;
    Expl.mwFreq=Expl.mwFreq*1e3; %GHz to MHz
    Expl.CrystalOrientation = [phi(iOrient) theta(iOrient) 0];
    % calculate energy levels and probabilites between them
    [EnergyLevels, Probabilities]=SolveEigenProblem(Hzf,zeemanop,Expl,Opt);
    
    % get populations
    %%% hey nino %%% - make user-supplied starting states possible
    Populations=diag(  sigeq(diag(EnergyLevels),Expl.Temperature)  );
    
    % transition frequencies between levels in MHz
    TransFreqs = abs(EnergyLevels-EnergyLevels');
    TransOffsets = (TransFreqs-Expl.mwFreq);
    
    % Calculate W1(HTA) at EDNMR position using experimental loaded Q-value
    GaussianPump=exp(-2*((TransFreqs-Expl.mwFreq)/Expl.ExciteWidth).^2);
%     GaussianPump=GaussianPump/max(GaussianPump(:));
    w1 = 2*pi*Expl.nu1.*GaussianPump;
    PumpWeights=(1-Bloch(w1,Expl.tHTA,Probabilities,Expl.Tm))/2; %between 0, and 1 (inversion)
    PumpWeights(PumpWeights<Opt.Threshold.Pump)=0;  % apply cutoff
    
    if ~any(PumpWeights) %abort this orientation if now transition is observed
        continue;
    end
    
    %calculate probe weight of each transition -  ensures orientation
    %selection
    ProbeWeights = Probabilities;
    ProbeWeights(ProbeWeights<Opt.Threshold.Probe)=0; %apply cutoff
        
    
    [Amp, Pos]=EDNMR_PeakList(TransFreqs,ProbeWeights,PumpWeights,Populations);
    
    %add up the spectrum from this powder orientation
    if ~isempty(Pos)
        %get rid of peaks slightly outside range
        %cannot really be done earlier, because it needs the knowledge of connected transitions)
        Amp(Pos<Expl.Range(1) | Pos>Expl.Range(2))=[];
        Pos(Pos<Expl.Range(1) | Pos>Expl.Range(2))=[];
        spec=spec+PowderWeight(iOrient)*makespec(Expl.Range,Expl.nPoints,Pos,Amp);
    end
    
end %end of powder loop


%build final spectrum
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

if nargin==3
    Opt=varargin{3};
else
    Opt=struct;
end
if ~isfield(Opt,'nKnots') Opt.nKnots=31; end
if ~isfield(Opt,'Threshold')
    Opt.Threshold=struct;
end
if ~isfield(Opt.Threshold,'Probe')  Opt.Threshold.Probe=1e-4; end
if ~isfield(Opt.Threshold,'Pump')  Opt.Threshold.Pump=1e-4; end
if ~isfield(Opt.Threshold,'Iso')  Opt.Threshold.Iso=1e-3; end
if ~isfield(Opt,'Output')
    Opt.Output='summed';
else
    if ~strcmp(Opt.Output,'summed') && ~strcmp(Opt.Output,'separate')
        error('Please specify Opt.Output as summed or separate ')
    end
end

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

%get cosine of flip anlge from weight
CosPump=-2*PumpWeights+1;

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
                        
                        Pos=[Pos TransFreqs(i_probe,j_probe)-TransFreqs(i_pump,j_pump)]; %postion of pumped transition wrt to observed one(not! the resonator. is this correct?)
                        
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

end %end of peak function