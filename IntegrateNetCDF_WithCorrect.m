function [RT,TC,X,IDrt,IDtc,IDx,IDrun,iX,xNat] = IntegrateNetCDF_WithCorrect(d,L,c)
% d is a string which is the file path to the .cdf (netcdf) file that will be integrated
% L is a string that is the name of the library (without the .m) used to integrate the data. The library can be in the MatLab directory, or in the folder containing the .cdf files to be integrated
% L can also be a function name (no single quotes, w/o .m) that is in the current MatLab director
% c is one if the data is to be corrected for natural abundances, zero if it is not to be corrected. 
% If c is one, startup.m must be ran before this function is run
% creates a .mat file and a .csv file with the outputs - the .mat file is
% call to function: IntegrateNetCDF_WithCorrect(d,L,c)

% Author: Maciek Antoniewicz (maciek@mit.edu)
% Date of last revision: 4/28/06 -- Maciek
% Also modified by Jamey Young?
% Also modified by Nathaniel Vacanti

iscorr=c;
directory=d;
MSLib=L;
%If the input MSLib is a string
%Include the directory of the data in the MatLab search path and
%Obtain the output of the Library file
if ischar(MSLib) == 1
    org_path=path; % save the original search path
    path(org_path,directory); % add the directory of the data to the search path
    Library = str2func(MSLib); % make the library name string a function with that name
    MSLib = Library(); % have the output of that function be the variable MSLib - MSLib is no longer a string
    path(org_path); % change the MatLab search path back to what it was
end

if nargin < 3
    iscorr = false;
end

% Init
RT = [];    % retention times of identified metabolites, and
IDrt = {};  % ... corresponding name of metabolite.
TC = [];    % integrated ion counts of selected fragments, and
IDtc = {};  % ... corresponding name of fragment.
X = [];     % relative intensities of selected ions, and
IDx = {};   % ... corresponding name of ion.
IDrun = {}; % name of the source file.
xNat = [];  % natural abundance distribution of fragments.

% Collect names of text-files in current folder
files = dir([directory,'\*.cdf']);

% Put file names in alphabetical order - added by Nate Vacanti 8/12
% put file names in cell array of strings
NumFiles = length(files);
FileNamesHolder = cell(NumFiles,1);
for i=1:NumFiles
    FileNamesHolder{i} = files(i).name;
end
FileNames = alphabetize(FileNamesHolder);

% Integrate GCMS data in each translated text-file
% instances of files(i).name replaced with FileNames{i} so integrate will
% work taking files in alphabetical order
for i=1:length(files)
    disp(['Analyzing data file: ',FileNames{i},' ...']);
    filename = FileNames{i};
    filename = filename(1:end-4);
    [rt,tc,x,IDrt,IDtc,IDx,xNat] = IntegrateGCMS_File([directory,'\',FileNames{i}],MSLib,i,['GC/MS data from file: ',filename],iscorr);
    RT = [RT,rt];
    TC = [TC,tc];
    X = [X,x];
    IDrun{end+1} = filename;
end

% Create list of indices for mass isotopomers of each fragment
iX = getIndex(MSLib);

% Write data to text file, that can then be imported into Excel (comma-separated)
Write2File(directory,RT,TC,X,IDrt,IDtc,IDx,IDrun,xNat,iscorr);
% Save results to Matlab-file
if iscorr == 0
    save([directory,'\Data(uncorrected).mat']);
elseif iscorr == 1
    save([directory,'\Data(corrected).mat']);
end


% -------------------------------------------------------------------------
function [rt,tc,x,IDrt,IDtc,IDx,xNat] = IntegrateGCMS_File(file,MSLib,figNr,figTitle,iscorr,tracer)

% Init
rt = [];    % retention times of identified metabolites, and
IDrt = {};  % ... corresponding name of metabolite
tc = [];    % integrated ion counts of selected fragments, and
IDtc = {};  % ... corresponding name of fragment
x = [];     % relative intensities of selected ions, and
IDx = {};   % ... corresponding name of ion
xNat = [];  % natural abundance distribution of selected fragments

% T = timepoints (in minutes)
% MS = ion counts for m/z 0-500 (rows - timepoints; columns - mass)
[T,MS] = readGCMSDataFile(file);
TIC = sum(MS');

% Integrate each fragment ion individually
for i=1:length(MSLib)
    
    % Locate the metabolite peak, and report the retention time
    [timePeak,idxPeak] = findTargetMetabolite(MSLib(i),T,MS);
    
    if timePeak % metabolite was found
        % Select area around the metabolite peak, 1 minute before and 1 minute after the peak
        idxArea = getSurroundingArea(idxPeak,T,1);
        idxPeak = find(idxArea==idxPeak);
        % Find the beginning and end of metabolite peak
        [peak,cnts,base,noise] = adjustPeakArea(idxPeak,MS(idxArea,:),MSLib(i).MainIons);
        % Integrate ion counts over the peak area
        [ix,itc] = integrateSIM(peak,MS(idxArea,:),MSLib(i).SelectedIons,base);
        [iIDx,iIDtc] = getIonNames(MSLib(i));
        
    else % metabolite was not found
        disp(['The following metabolite was not found: ',MSLib(i).Name])
        [ix,itc] = integrateSIM([],[],MSLib(i).SelectedIons,[]);
        [iIDx,iIDtc] = getIonNames(MSLib(i));
        
    end
    
    % Calculate theoretical natural abundances
    xNat = [xNat ; naturalSIM(MSLib(i).SelectedIons,MSLib(i).SelectedIonsFormula)];
    
    % Collect data in vectors
    
    if iscorr
        allatoms = txt2atm(MSLib(i).SelectedIonsFormula);
        %labatoms = txt2atm(MSLib(i).SelectedIonsAtoms);
        %more = allatoms-labatoms;
        ix = mscorrect(ix,allatoms); %,labatoms);
    end
    
    x = [x;ix];
    IDx = [IDx;iIDx];
    tc = [tc;itc];
    IDtc = [IDtc;iIDtc];
    rt = [rt;timePeak];
    IDrt{end+1,1} = MSLib(i).Name;
    
end


% -------------------------------------------------------------------------
function [t_rd matrix] = readGCMSDataFile(file)
%Written by Nathaniel Vacanti
%special thanks to Ryan LaCroix

% Get times, masses, intensities, and index for time separation from cdf file
t = ncread(file,'scan_acquisition_time')/60;
t_rd=round(100*t)/100;
t_rd=t_rd';
m = ncread(file,'mass_values');
int = ncread(file,'intensity_values');
ind = ncread(file,'scan_index');

m_rd_pre = round(m*100)/100; %first round to the nearest 0.01
m_rd = round(m_rd_pre-0.15); %now round to nearest whole number
%rounding twice done to mimmick previous method of writing .cdf file to
%.txt file first - can be omitted


% create a column with 500 integer masses from one to 500, if reading from
% amu 1 to 500 or less, else read to highest mass
if max(m_rd)<=5000;
    m_srt = 1:1:5000;
else
    m_srt = 1:1:max(m_rd);
end
m_srt = m_srt';

matrix = zeros(length(t), length(m_srt));


for i=1:length(t) %iterate over time
    if i==length(t) %if you are at the last time point
        v=length(m_rd); %read to the end of the rounded mass vector
    else
        v=ind(i+1); %otherwise read until the next index point
    end
    for k=ind(i)+1:v %read the mass vector from the first point of the new index, to the next index point
        matrix(i,m_rd(k))=int(k); %fill in matrix with intensities corresponding to masses at that time
        %The column number is the mass - only true because mass starts at
        %one and goes by one to end. Very convenient, saves me from having
        %to search through mass vector 500+ (the # of mass points) times
    end
    
end

% -------------------------------------------------------------------------
function [timepoint,idx] = findTargetMetabolite(MSLibItem,T,MS)

% Author: Maciek Antoniewicz (maciek@mit.edu)
% Date of last revision: 3/26/06 -- Maciek

% Collect info about target metabolite from library
RetTime = MSLibItem.RetTime;
MainIons = MSLibItem.MainIons;
MainIonsRelInt = MSLibItem.MainIonsRelInt;
    
% Size of data matrix
[nT,nM] = size(MS);

% Rescale relative intensities
nMainIons = size(MainIons,1);
MainIonsRelInt = MainIonsRelInt/sum(MainIonsRelInt);
A = zeros(nMainIons,nM);
for i=1:nMainIons
    A(i,MainIons(i,1):MainIons(i,2)) = 1;
end
sim = A*MS';
relsim = sim./(repmat(sum(sim),nMainIons,1)+1);

% Determine most likely location of peak based on relative peak values
K = nMainIons*0.1;
p1 = (sum(abs(relsim-repmat(MainIonsRelInt,1,nT)))/K+1).^(-2);

% Determine most likely location of peak based on retention time
K = 0.05; % minutes
p2 = normpdf(T-RetTime,0,K);

% Determine most likely location of metabolite peak
if length(p1)==1
    x = sim.*p2;
else
    x = sum(sim).*p1.*p2;
end
idx = find(x==max(x));
idx = idx(1);
timepoint = T(idx);

% Set a minimum threshold for metabolite identification
if size(MainIons,1) > 1
    if p1(idx)<0.1 | p2(idx)<1e-4 % metabolite was not found
        timepoint = 0;
        idx = [];
    end
end


% -------------------------------------------------------------------------
function idxArea = getSurroundingArea(idxPeak,T,dT);

% Determine scan frequency in the area around the peak 
idx = max(1,idxPeak-10):min(length(T),idxPeak+10);
ScansPerMin = (length(idx)-1)/(T(idx(end))-T(idx(1)));

% Select area around the peak
NScans = round(dT*ScansPerMin);
idxArea = max(1,idxPeak-NScans):min(length(T),idxPeak+NScans);
idxPeak = find(idxArea==idxPeak);


% -------------------------------------------------------------------------
function [peak,cnts,base,noise] = adjustPeakArea(idxPeak,MS,SelectedIons)
% This function finds the beginning and end of a metabolite peak.
% 
% INPUT:
% idxPeak = scan number at which the metabolite peak was identified
% MS = matrix with ion counts (rows = scan number; columns = mass 0-500)
% SelectedIons = the selected mass range(s) for integration later on
% SelectedIons = [mass_range_begin_1,mass_range_begin_1 ; 
%                 mass_range_begin_2,mass_range_begin_2 ; ...]  
%
% OUTPUT:
% peak = scan numbers of beginning and end of metabolite peak
% cnts = total ion counts of the selected mass range(s)
% base = scan numbers that were identified as being at baseline level
% noise = average intensity of noise in the peak area

% Author: Maciek Antoniewicz (maciek@mit.edu)
% Date of last revision: 3/26/06 -- Maciek

% Determine total ion counts for the selected mass range(s)
mz = [];
for i=1:size(SelectedIons,1)
    mz = [mz,SelectedIons(i,1):SelectedIons(i,2)];
end
if length(mz)==1 
    cnts = MS(:,mz)';
else 
    cnts = sum(MS(:,mz)');
end

% Smooth out data
[c,dc] = QuadSlope(1:length(cnts),cnts,2);

% Estimate average noise-intensity (minimum noise level is also set);
NoiseLB = 10;
noise = max(NoiseLB,1.8*median(abs(c-cnts)));
    
% Find scans where the intensity is at baseline level
idxBaseline = findBaselineAreas(c,noise,4);

% Find average baseline level nearest to the peak
base = findNearestBaseline(idxPeak,idxBaseline);
bsln = mean(cnts(base));

% Init
% Initial scan to determine edges of the peak
THinit = 0.3;   % Include all adjesent scans that are >30% of peak height
LB = max(1,idxPeak-1);
while LB>1 & cnts(LB)>THinit*cnts(idxPeak)
    LB = LB-1;
end
UB = min(length(cnts),idxPeak+1);
while UB<length(cnts) & cnts(UB)>THinit*cnts(idxPeak)
    UB = UB+1;
end
peak = [LB,UB];

% STEP 1. First adjust the beginning of the peak
% Thresholds for determination when the peak ends
THc = 10; % how many times the intensity needs to be higher than noise level
THdc = 2; % how many times the increase in intensity from scan-to-scan needs to be higher than noise level
ready = 0;
idx = peak(1);
if dc(idx)>THdc*noise & c(idx)-bsln>THc*noise
    % Decrease lower bound
    while ~ready
        if idx==1
            ready = 1;
        else
            idx = idx-1;
        end
        if dc(idx)<THdc*noise | c(idx)-bsln<THc*noise
            ready = 1;
            peak(1) = idx+1;
        end
    end
elseif dc(idx)<THdc*noise & c(idx)-bsln>THc*noise
    % Increase lower bound
    while ~ready
        if idx==peak(2)
            ready = 1;
        else
            idx = idx+1;
        end
        if dc(idx)>THdc*noise | c(idx)-bsln<THc*noise
            ready = 1;
            peak(1) = idx;
        end
    end
else
    ready = 1;
end

% STEP 2. Next adjust the end of the peak
ready = 0;
idx = peak(2);
if dc(idx)<-THdc*noise & c(idx)-bsln>THc*noise
    % Increase upper bound
    while ~ready
        if idx==length(c)
            ready = 1;
        else
            idx = idx+1;
        end
        if dc(idx)>-THdc*noise | c(idx)-bsln<THc*noise
            ready = 1;
            peak(2) = idx-1;
        end
    end
elseif dc(idx)>-THdc*noise & c(idx)-bsln>THc*noise
    % Decrease upper bound
    while ~ready
        if idx==peak(1)
            ready = 1;
        else
            idx = idx-1;
        end
        if dc(idx)<-THdc*noise | c(idx)-bsln<THc*noise
            ready = 1;
            peak(2) = idx;
        end
    end
else
    ready = 1;
end


% -----------------------------------------------------------------------
function [yfit,dydx] = QuadSlope(x,y,varargin)

n = 2; % default is quadratic approximation
if nargin > 2
    n = varargin{1};
end

for i=1:length(y)-2*n
    p = FitSlope(x(i:i+2*n),y(i:i+2*n),n);
    if i==1
        for j=1:n
            yfit(j) = p(1)*x(j)^2 + p(2)*x(j) + p(3);
            dydx(j) = 2*p(1)*x(j) + p(2);
        end
    elseif i==length(y)-2*n
        for j=i+(n+1:2*n)
            yfit(j) = p(1)*x(j)^2 + p(2)*x(j) + p(3);
            dydx(j) = 2*p(1)*x(j) + p(2);
        end
    end
    yfit(i+n) = p(1)*x(i+n)^2 + p(2)*x(i+n) + p(3);
    dydx(i+n) = 2*p(1)*x(i+n) + p(2);
    dydx2(i+n) = 2*p(1); 
end


% -----------------------------------------------------------------------
function p = FitSlope(x,y,n)

A = zeros(3);
b = zeros(3,1);
for i=1:2*n+1
    A(1,1) = A(1,1) + x(i)^4;
    A(1,2) = A(1,2) + x(i)^3;
    A(1,3) = A(1,3) + x(i)^2;
    A(2,1) = A(2,1) + x(i)^3;
    A(2,2) = A(2,2) + x(i)^2;
    A(2,3) = A(2,3) + x(i);
    A(3,1) = A(3,1) + x(i)^2;
    A(3,2) = A(3,2) + x(i);
    A(3,3) = A(3,3) + 1;
    b(1,1) = b(1,1) + x(i)^2*y(i);
    b(2,1) = b(2,1) + x(i)*y(i);
    b(3,1) = b(3,1) + y(i);
end
p = pinv(A)*b;


% -----------------------------------------------------------------------
function idxBaseline = findBaselineAreas(data,noise,n);

idx = [];
while isempty(idx)
    % Find potential baseline locations
    idx1 = find(abs(data(1:end-1)-data(2:end))<=noise);
    idx2 = find(abs(data(1:end-2)-data(3:end))<=noise);
    idx3 = find(abs(data(1:end-3)-data(4:end))<=noise);
    idx4 = find(abs(data(1:end-4)-data(5:end))<=noise);
    idx5 = find(abs(data(1:end-5)-data(6:end))<=noise);
    idx =  intersect(intersect(intersect(intersect(idx1,idx2),idx3),idx4),idx5);
    % Adjust assumed noise level to find a baseline    
    if isempty(idx)
        noise = 2*noise;
    end
end

% Init
baseline = [];

if length(idx)

    % Fill in single gaps between baseline locations
    idx = sort([idx,idx(find(idx(2:end)-idx(1:end-1)==2))+1]);
    endIdx = find(idx(2:end)-idx(1:end-1)~=1)';
    if isempty(endIdx)
        endIdx = length(idx);
    elseif endIdx(end) ~= length(idx)
        endIdx(end+1) = length(idx);
    end
    beginIdx = [1;endIdx(1:end-1)+1];
    
    % Make column vectors
    beginIdx = beginIdx(:);
    endIdx = endIdx(:);
    
    % Locations of baseline
    baseline = [idx(beginIdx)',idx(endIdx)'+5];

    % Select only those baseline that are at least 'n' timepoints long
    if any(endIdx-beginIdx+1>=n)
        baseline = baseline(find((endIdx-beginIdx+1)>=n),:);
    else
        maxLengthBaseline = max(endIdx-beginIdx+1);
        baseline = baseline(find((endIdx-beginIdx+1)>=maxLengthBaseline),:);
    end
    
end

% Index to baseline timepoints
idxBaseline = [];
for i=1:size(baseline,1)
    idxBaseline = [idxBaseline,baseline(i,1):baseline(i,2)];
end


% -----------------------------------------------------------------------
function base = findNearestBaseline(idxPeak,idxBase)

% Find nearest 40 scans of baseline that are closest to the peak
base = [abs(idxBase-idxPeak);1:length(idxBase)];
base = sortrows(base')';
if size(base,2)>40
    base = base(2,1:40);
else
    base = base(2,:);
end
base = idxBase(base);


% -------------------------------------------------------------------------
function [x,tc] = integrateSIM(peak,MS,SIM,base)

% Author: Maciek Antoniewicz (maciek@mit.edu)
% Date of last revision: 3/26/06 -- Maciek

% Init
x = [];     % Relative ion abundances
tc = [];    % Total counts per fragment

if ~isempty(peak) & (peak(1)~=peak(2))
    % Integrate SIM data
    sumMS = sum(MS(peak(1):peak(2),:));
    for i=1:size(SIM,1)
        % Calculate total counts
        simCnts = sumMS(SIM(i,1):SIM(i,2));
        % Calculate baseline counts
        baseCnts_median = median(MS(base,SIM(i,1):SIM(i,2)))*(peak(2)-peak(1)+1);
        baseCnts_mean   = mean(MS(base,SIM(i,1):SIM(i,2)))*(peak(2)-peak(1)+1);
        baseCnts = max([baseCnts_mean;baseCnts_median]);
        % Correct for baseline
        simCntsCorrected = max(0,simCnts-baseCnts);
        % Calculate total integrated counts
        tc = [tc;sum(simCntsCorrected)];
        % Convert counts into abundances
        if sum(simCntsCorrected)
            if 1 % Convert counts into relative abundances
                x = [x;simCntsCorrected'/sum(simCntsCorrected)];
            else % Convert counts into abundances scaled to largest abundance
                x = [x;100*simCntsCorrected'/max(simCntsCorrected)];
            end
        else
            x = [x;simCntsCorrected'];
        end
        
    end
    
else % metabolite was not found
    for i=1:size(SIM,1)
        x = [x;zeros(SIM(i,2)-SIM(i,1)+1,1)];
        tc = [tc;0];
    end
end


% -------------------------------------------------------------------------
function [MIDV,ExactMass,Composition] = MSNatural(compound,varargin)

% INPUT: 
%   compound        - either a string with chemical formula (e.g. 'C12H34O5N2')
%                   - or a vector (size 1x32) with elemental composition:
%                       v(1)  = number of H-atoms
%                       v(12) = number of C-atoms
%                       v(14) = number of N-atoms
%                       v(16) = number of O-atoms
%                       v(28) = number of Si-atoms
%                       v(28) = number of Si-atoms
%                       v(31) = number of P-atoms
%                       v(32) = number of S-atoms
%
%   varargin{1}     - preferred length of simulated MIDV (default n=15).
%
% OUTPUT: 
%   MIDV            - mass isotopomer distribution vector with the natural abundance of given compound
%   ExactMass       - exact mass of the M+0 ion of given compound

% Author: Maciek Antoniewicz
% Date: 07/12/04

% Database of natural isotopes
DB = IsotopeDatabase;

% Determine type of input (string or vector)
Composition = zeros(1,32);
if isa(compound,'double')
    if size(compound,1)==1 & size(compound,2)==32
        Composition = compound;    
    else
        error('Error: invalid input, vector should have length 32.')
    end
elseif isa(compound,'char')
    Composition = ReadComposition(compound,DB);
else
    error('Error: invalid input')
end

% Prepare index vectors used for construction of matrix A
n = 15; 
if nargin > 1
    n = varargin{1};    
end
v = (1:n)';
V = v(:,ones(1,n));
col = V';
row = col + [zeros(1,n);V(1:n-1,:)];

% Init
ExactMass = 0;
MIDV = zeros(n,1);
MIDV(1) = 1;

% Construct isotope distribution and calculate exact mass of compound
for i=1:length(Composition)
    if Composition(i)
        % Construct matrix A
        NAtoms = Composition(i);
        Y = DB(i).NatAbund';
        if length(Y) > n
            Y = Y(1:n);
        end
        Ly = length(Y);
        A = sparse(row(1:Ly,1:n),col(1:Ly,1:n),Y(:,ones(n,1)));
        A = A(1:n,1:n);
        % Update isotope distribution and exact mass
        MIDV = (A^NAtoms)*MIDV;
        ExactMass = ExactMass + DB(i).ExactMass(1)*NAtoms;
    end
end
    

% ---------------------------------------------------------------------------------------
function Composition = ReadComposition(Name,DB)

% Numbers = 48 - 57
% Space = 32
% Small letters = 97 - 122
% Capital letters = 65 - 90

i = 1;
Composition = zeros(1,32);
while i<=length(Name)
    S1 = '';
    S2 = '';
    
    % Read first field
    if double(Name(i))==32
    
        % Skip field
        i = i+1;
        
    else
        
        if double(Name(i))>=97 & double(Name(i))<=122
            S1 = char(double(Name(i))-32);
        elseif double(Name(i))>=65 & double(Name(i))<=90
            S1 = Name(i);
        elseif double(Name(i))>=48 & double(Name(i))<=57
            error('Error: reading metabolite composition')
        end
        
        % Read second field
        if i<=length(Name)-1
            if double(Name(i+1))>=97 & double(Name(i+1))<=122
                S2 = [S1,Name(i+1)];
            elseif double(Name(i+1))>=65 & double(Name(i+1))<=90
                S2 = [S1,char(double(Name(i+1))+32)];
            elseif double(Name(i+1))>=48 & double(Name(i+1))<=57
                % do nothing
            end
        end
    
        % Determine element. First try long name.
        Element = find(strcmp({DB.Element},S2));
        if isempty(Element)
            Element = find(strcmp({DB.Element},S1));
            if isempty(Element)
                error('Error: reading metabolite composition')
            else
                i=i+1;
            end
        else
            i=i+2;
        end
    
        % Determine number of atoms
        if i<=length(Name)
            stop = 0;
            NumAtoms = '';
            while ~stop & i<=length(Name)
                if double(Name(i))>=48 & double(Name(i))<=57
                    NumAtoms = [NumAtoms,Name(i)];
                    i=i+1;
                else
                    stop = 1;
                end
            end
            if ~isempty(NumAtoms)
                NumAtoms = str2num(NumAtoms);    
            else
                NumAtoms = 1;
            end
        else
            NumAtoms = 1;
        end
    
        % Update composition
        Composition(Element) = Composition(Element) + NumAtoms;
    
    end
    
end


% ---------------------------------------------------------------------------------------
function DB = IsotopeDatabase()

DB(1).Element = 'H';
DB(1).ExactMass = [1.007825,2.014102];
DB(1).NatAbund = [99.9844,0.0156]/100;
DB(1).Z = 1;

DB(12).Element = 'C';
DB(12).ExactMass = [12,13.003355];
DB(12).NatAbund = [98.9184,1.0816]/100; % Corresponding to d13C vs. VPDB = -27 o/oo
DB(12).Z = 6;

DB(14).Element = 'N';
DB(14).ExactMass = [14.003074,15.000109];
DB(14).NatAbund = [99.6337,0.3663]/100;
DB(14).Z = 7;

DB(16).Element = 'O';
DB(16).ExactMass = [15.994915,16.999132,17.999160];
DB(16).NatAbund = [99.758,0.038,0.204]/100;
DB(16).Z = 8;

DB(28).Element = 'Si';
DB(28).ExactMass = [27.976927,28.976495,29.973770];
DB(28).NatAbund = [92.22,4.69,3.09]/100;
DB(28).Z = 14;

DB(31).Element = 'P';
DB(31).ExactMass = [30.973762];
DB(31).NatAbund = [100]/100;
DB(31).Z = 15;

DB(32).Element = 'S';
DB(32).ExactMass = [31.972071,32.971458,33.967867,35,35.967081];
DB(32).NatAbund = [95.039,0.749,4.197,0,0.015]/100;
DB(32).Z = 16;


% -------------------------------------------------------------------------
function [IDx,IDtc] = getIonNames(MSLibItem)

% Init
IDx = {};   % Ion names
IDtc = {};  % Fragment names

% Collect info about metabolite
Metabolite = MSLibItem.Name;
SIM = MSLibItem.SelectedIons;

% Create list with fragment and ion names
for i=1:size(SIM,1)
    ions = SIM(i,1):SIM(i,2);
    IDtc{end+1,1} = [Metabolite,' ',int2str(ions(1))];
    for j=1:length(ions)
        IDx{end+1,1} = [Metabolite,int2str(ions(j)),' (M',int2str(j-1),')'];
    end
end


% -------------------------------------------------------------------------
function data = naturalSIM(SIM,Formula)

% Calculate natural isotope abundance distribution
data = [];
for i=1:size(SIM,1)
    natural = MSNatural(Formula{i});
    Lsim = SIM(i,2)-SIM(i,1)+1;
    if length(natural) >= Lsim
        % Return relative abundances
        if 1
            data = [data;natural(1:SIM(i,2)-SIM(i,1)+1)];
        % Return abundances scaled to largest abundance
        else
            data = [data;100*natural(1:SIM(i,2)-SIM(i,1)+1)/max(natural(1:SIM(i,2)-SIM(i,1)+1))];
        end
        
    else
        natural = [natural;zeros(Lsim,1)];
        % Return relative abundances
        if 1
            data = [data;natural(1:SIM(i,2)-SIM(i,1)+1)];
        % Return abundances scaled to largest abundance
        else
            data = [data;100*natural(1:SIM(i,2)-SIM(i,1)+1)/max(natural(1:SIM(i,2)-SIM(i,1)+1))];
        end
        
    end
end


% -------------------------------------------------------------------------
function iX = getIndex(MSLib);

% Create list of indices for mass isotopomers of each fragment
iX = {};
cnt = 0;
for i=1:length(MSLib)
    for j=1:size(MSLib(i).SelectedIons,1)
        L = MSLib(i).SelectedIons(j,2) - MSLib(i).SelectedIons(j,1) + 1;
        iX{end+1} = [cnt+1:cnt+L];
        cnt = cnt + L;
    end
end


% -------------------------------------------------------------------------
function outputfile = Write2File(directory,RT,TC,X,IDrt,IDtc,IDx,IDrun,xNat,iscorr)

% Create output file
if iscorr==0
    outputfile = [directory,'\Data(uncorrected).csv'];
elseif iscorr==1
    outputfile = [directory,'\Data(corrected).csv'];
end
fid = fopen(outputfile,'w');

% Create headers
fprintf(fid,'Retention Times (min.)\n');
fprintf(fid,'Metabolite');
for i=1:length(IDrun)
    fprintf(fid,[',',IDrun{i}]);
end
fprintf(fid,'\n');
% Write retention times
for i=1:size(RT,1)
    fprintf(fid,[IDrt{i},',']);
    for j=1:size(RT,2)
        fprintf(fid,'%2.2f',RT(i,j));
        fprintf(fid,',');            
    end
    fprintf(fid,'\n');
end

% Create headers
fprintf(fid,'\nTotial ion counts\n');
fprintf(fid,'Fragment');
for i=1:length(IDrun)
    fprintf(fid,[',',IDrun{i}]);
end
fprintf(fid,'\n');
% Write total ion counts
for i=1:size(TC,1)
    fprintf(fid,[IDtc{i},',']);
    for j=1:size(TC,2)
        fprintf(fid,'%6.0f',TC(i,j));
        fprintf(fid,',');          
    end
    fprintf(fid,'\n');
end

% Create headers
fprintf(fid,'\nMass isotopomer abundances\n');
fprintf(fid,'Ion');
for i=1:length(IDrun)
    fprintf(fid,[',',IDrun{i}]);
end
fprintf(fid,',Theory, \n');
% Write relative abundances
for i=1:size(X,1)
    if i > 1 && isequal(IDx{i}(end-3:end),'(M0)')
        fprintf(fid,'\n');
    end
    fprintf(fid,[IDx{i},',']);
    for j=1:size(X,2)
        fprintf(fid,'%0.4f',X(i,j));
        fprintf(fid,',');            
    end
    fprintf(fid,'%0.4f',xNat(i));
    fprintf(fid,'\n');
end

% Finish writing to output file
fclose(fid);

%--------------------------------------------------------------------------

function [x,M] = mscorrect(x,y)

%MSCORRECT Correct mass isotopomer distribution for natural abundance.
%   X = MSCORRECT(X,Y) corrects the MIDV in X for naturual isotopic
%   abundances. The input Y is a string or vector that gives the fragment
%   composition. 
%
%   X = MSCORRECT(X,Y,Z) assumes that any deviation from natural abundance
%   is due to incorporation of labeled Z atoms. This is desirable if Z is a
%   tracer atom. Uses method of Fernandez et al. (1996) Journal of Mass
%   Spectrometry 31 pp. 255--262.
%
%   [X,M] = ... also returns the correction matrix.

% determine tracer identity
% Match labeled atom (in tracer) to number below in y(?)
% z=1 for H
% z=6 for C
% z=7 for N
% z=8 for O
% z=9 for F
% z=14 for Si
% z=15 for P
% z=16 for S
z = 6;

y = txt2atm(y);
ncarbons = y(z);
lorig = length(x);
l = min(length(x),ncarbons+1);
x = x(1:l);
M = zeros(l,l);
for i = 1:size(M,1)
    M(i,:) = [zeros(1,i-1),msnatural(y,l-i+1)'];
    y(z) = y(z)-1;
    if y(z) < 0
        break
    end
end
x = (x'*inv(M))';
x = [x;zeros(max(0,lorig-l),1)];

%--------------------------------------------------------------------------

function Alphabetized = alphabetize(ArrayOfStrings)
% Input - cell array of strings
% Output - the cell array of strings alphabetized
   % Note the Matlab sort function treats all capital letters as coming
   % before all lowercase letters - this function does not
   % returns a row cell array or column cell array depending on the input
   % cannot take a matrix cell array as an input


NumStrings = length(ArrayOfStrings);
ArrayLower = lower(ArrayOfStrings); %you need to do this before sorting because matlabe arranges capital letters before lowercase
[ArrayLower ix] = sort(ArrayLower);

%Initialize Alphabetized
[r c] = size(ArrayOfStrings);

if c==1
    Alphabetized = cell(NumStrings,1);
elseif r==1
    Alphabetized = cell(1,NumStrings);
end

for i=1:NumStrings % this loop recovers capital letters
    Alphabetized{i}=ArrayOfStrings{ix(i)};
end

%--------------------------------------------------------------------------

function y = txt2atm(s)
%TXT2ATM Convert a string representation of a chemical formula into
%its associated vector format.
%   Y = TXT2ATM(S) converts the string S into its composition vector
%   Y. S consists of a sequence of chemical symbols and numbers such 
%   as C6H12O6. Y(I) contains the number of atoms of the element with
%   atomic number I. By default, Y has length 16.

% set default parameters
ny = 16; % default length of composition vector

% initialize output
y = zeros(ny,1);

% check input
if isnumeric(s), % if the input is already a vector
    y(1:length(s)) = s; % return its current value
    return
elseif ischar(s) || iscell(s),
    s = cellstr(s);
    s = sprintf('%s',s{:});
else % unexpected data type
    error('Invalid input.')
end

% parse through the input string, one element at a time
while ~isempty(s), % loop while the string is not empty
    [s,z] = element(s); % determine element name
    [s,n] = getnumber(s); % determine number of atoms
    y(z) = y(z) + n; % update composition
end

% --------------------------------------------------------------------
function [s,z] = element(s)

z = []; % init
s = strtrim(s); % skip blanks
if isempty(s), return, end

% determine element, first try long name
z = dbquery(s(1:min(2,end)),'z');
if isempty(z), 
    z = dbquery(s(1),'z');
    if isempty(z), 
        error('Cannot read metabolite composition.')
    else
        s(1) = [];
    end
else 
    s(1:min(2,end)) = [];
end

% --------------------------------------------------------------------
function [s,n] = getnumber(s)

n = 1; % init
s = strtrim(s); % skip blanks
if isempty(s), return, end

% read in numbers until a non-number is found
N = '';
while ~isempty(s) && ismember(s(1),'0123456789'),
    N(end+1) = s(1);
    s(1) = [];
end
if ~isempty(N), n = sscanf(N,'%u'); end

%--------------------------------------------------------------------------

function varargout = dbquery(s,varargin)
%DBQUERY Query database of stable isotopes.
%   [V1,V2,...] = DBQUERY(S,P1,P2,...) returns a vector that contains
%   information on the stable isotopes of S. S can be a string containing
%   the atomic symbol or the atomic number of the element. The outputs
%   V1, V2, etc. correspond to the properties P1, P2, etc. where valid
%   properties are "s" (element symbol), "m" (vector of exact masses),
%   "n" (vector natural abundances) or "z" (atomic number).

% load database
persistent x
if isempty(x),
    x(1).s = 'H'; % atomic symbol
    x(1).m = [1.007825;2.014102]; % exact mass of each isotope
    x(1).n = [99.9844;0.0156]/100; % natural abundance of each isotope
    x(1).z = 1; % atomic number

    x(6).s = 'C';
    x(6).m = [12;13.003355];
    x(6).n = [98.9184;1.0816]/100; 
    % corresponding to d13C vs. VPx = -22 o/oo
    x(6).z = 6;

    x(7).s = 'N';
    x(7).m = [14.003074;15.000109];
    x(7).n = [99.6337;0.3663]/100;
    x(7).z = 7;

    x(8).s = 'O';
    x(8).m = [15.994915;16.999132;17.999160];
    x(8).n = [99.758;0.038;0.204]/100;
    x(8).z = 8;

    x(9).s = 'F';
    x(9).m = [18.998403];
    x(9).n = [100]/100;
    x(9).z = 9;

    x(14).s = 'Si';
    x(14).m = [27.976927;28.976495;29.973770];
    x(14).n = [92.22;4.69;3.09]/100;
    x(14).z = 14;

    x(15).s = 'P';
    x(15).m = [30.973762];
    x(15).n = [100]/100;
    x(15).z = 15;

    x(16).s = 'S';
    x(16).m = [31.972071;32.971458;33.967867;35;35.967081];
    x(16).n = [95.039;0.749;4.197;0;0.015]/100;
    x(16).z = 16;
end

% handle string input
if ischar(s), s = strcmpi({x.s},s); end

% assign outputs
len = length(varargin);
varargout = cell(1,len);
if ~any(s), return, end
for i = 1:len,
    varargout{i} = x(s).(varargin{i});
end

%--------------------------------------------------------------------------

function [x,m] = msnatural(y,varargin)
%MSNATURAL Calculate natural isotopic labeling.
%   [X,M] = MSNATURAL(Y) finds the natural isotopic labeling of the
%   fragment Y. Y can be a string or composition vector. X is the MIDV of
%   the naturally labeled fragment and M is its base mass.
%
%   [X,M] = MSNATURAL(Y,N) only computes elements of X for M+0 to M+N-1.

% convert input, if needed
y = txt2atm(y);

% determine length of midv to return
if nargin > 1, 
    nx = varargin{1};
else
    nx = maxmz(y) + 1;
end

% prepare index vectors for construction of matrix A
V = tile((1:nx)',1,nx);
col = V';
row = col + [zeros(1,nx);V(1:nx-1,:)];

% init
m = 0;
x = zeros(nx,1);
x(1) = 1;

% construct mass isotopomer distribution and calculate exact mass of compound
for i = find(y)',
    
    % construct matrix A
    [M,Y] = dbquery(i,'m','n');
    if length(Y) > nx,
        Y = Y(1:nx);
    end
    nY = length(Y);
    A = sparse(row(1:nY,1:nx),col(1:nY,1:nx),Y(:,ones(nx,1)));
    A = A(1:nx,1:nx);
    
    % update isotope distribution and exact mass
    x = full((A^y(i))*x);
    m = m + M(1)*y(i);
    
end

% remove insignificant trailing zeros if desired length not specified
if nargin < 2, x = x(x >= eps); end 

%--------------------------------------------------------------------------

function B = tile(A,M,N)
%TILE Replicates a matrix.
%   B = TILE(A,M,N) forms B by stacking M-by-N copies of A. 

% Created by Jamey D. Young on 30-Nov-2005.

[m,n] = size(A); % get dimensions
nr = m*M; % store number of rows of B
nc = n*N; % store number of columns of B
i = 1 + rem(0:nr-1,m); % compute row indices
j = 1 + rem(0:nc-1,n); % compute column indices
B = A(i,j); % tile matrix
