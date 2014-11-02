%rainSTORM_2BPolarizationAnalysis
%This script uses the main RainSTORM Software to analyse 2 stacks of x-y
%polarized images

% Made by Daan van Kleef & Mehdi Goudarzi

% Only change to false if analysing multiple angles for plotting
flagInitialRun = true;

if flagInitialRun 
    clear
    flagInitialRun = true;
end

%Real angle of polarization
Phireal = 0;

% Flagsetzero sets missing observations to count 10^-10
% Flagsetbackground sets missing observations to background level
% Flagbackgroundsubtraction subtracts background level from observed photon
% counts and puts missing observations to 0
% NOTE, do not use backgroundsubtraction and setbackground together!
Flagsetzero = 0;
Flagsetbackground = 1;
Flagbackgroundsubtraction = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Use rainSTORM to analyse raw image data
%    Then read the data from the Matlab struct into a useful array...
startup  % Run rainSTORM once (manually), on the X-polarisation data
uiwait 

% Read in as an array!
Nxtable = struct2table((params.localization.results.SupResParams));
for count = 1:size(Nxtable,2)
    Nxdummy = double(table2array(Nxtable(:,count)));
    if count == 1
        Nxarray = Nxdummy;
    else
        Nxarray = cat(2,Nxarray,Nxdummy);
    end
end

startup  % Run rainSTORM once (manually), on the Y-polarisation data
uiwait

% Read in as an array!
Nytable = struct2table((params.localization.results.SupResParams));
for count = 1:size(Nytable,2)
    Nydummy = double(table2array(Nytable(:,count)));
    if count == 1
        Nyarray = Nydummy;
    else
        Nyarray = cat(2,Nyarray,Nydummy);
    end
end




% 2. Identify matching images
%    Method 1: assume fluorophore position is known (from simulation
%    params)

% Method 1:
% Remove spurious localizations from Nxarray and place in Nxarraycorrect
% Place rejected localizations into Nxreject
NxarrayCorrect = Nxarray;
NxDistance     = sqrt( (Nxarray(:,2)-32).^2 + (Nxarray(:,3)-30).^2 );
Nxreject = NxarrayCorrect(NxDistance > 1, :);
NxarrayCorrect(NxDistance > 1, :) = [];

%Idem for Nyarray
NyarrayCorrect = Nyarray;
NyDistance     = sqrt( (Nyarray(:,2)-32).^2 + (Nyarray(:,3)-30).^2 );
Nyreject = NyarrayCorrect(NyDistance > 1, :);
NyarrayCorrect(NyDistance > 1, :) = [];

numberOfFrames = params.localization.results.numberOfFrames;

%Initiate matrix to hold combined count list
% ComCount will hold the combined count list.
% Column Structure:
    %  Frame ID  | Nx | Ny | Nx present? | Ny present? | Estimated angle
    % Last 2 column will hold flags to indicate succesful counts
        % 0 = absent, 1 = present
        % Can be used to relocalize missing observations later
ComCount = zeros(numberOfFrames,6);

%Initiate frame ID list
for count = 1:numberOfFrames
    ComCount(count,1) = count;
end

%Copy data from Nxarraycorrect into ComCount
for count = 1:size(NxarrayCorrect)
    FrameID = NxarrayCorrect(count,1);
    ComCount(FrameID,2) = NxarrayCorrect(count, 5);
    ComCount(FrameID,4) = 1; 
end

%Copy data from Nyarraycorrect into Comcount
for count = 1:size(NyarrayCorrect)
    FrameID = NyarrayCorrect(count,1);
    ComCount(FrameID,3) = NyarrayCorrect(count, 5);
    ComCount(FrameID,5) = 1; 
end

%2b. (optional) Determine average residual background level from rejected
%localisations

ResBackground = 0.5*(mean(Nxreject(:,5)) + mean(Nyreject(:,5)));



% 3. Complete missing paired observations

% Method 1: Assume missing observations are count 10^-10
if(Flagsetzero == 1);

for count = 1:numberOfFrames
    if ComCount(count, 4) == 0
        ComCount (count , 2) = 10^(-10);
    end
    if ComCount(count,5) == 0
        ComCount(count,3) = 10^(-10);
    end
end
end

% Method 2: Assume missing observations are background count level
if(Flagsetbackground == 1);
for count = 1:numberOfFrames
    if ComCount(count, 4) == 0
        ComCount (count , 2) = ResBackground;
    end
    if ComCount(count,5) == 0
        ComCount(count,3) = ResBackground;
    end
end
end

% Method 3. Subtract average residual background level from measured
%accepted intensities and put missing observations to 0
if Flagbackgroundsubtraction == 1
    ComCount(:,2:3) = ComCount(:,2:3)-ResBackground;
    for count = 1:numberOfFrames
    if ComCount(count, 4) == 0
        ComCount (count , 2) = 10^(-10);
    end
    if ComCount(count,5) == 0
        ComCount(count,3) = 10^(-10);
    end
end
end
% 4. Caclulate polarization angle
ComCount(:,6) = (180/pi)*acot(sqrt(ComCount(:,2)./ComCount(:,3)));


%5. Conclusions
% Results is a vector with the avergae estimated angle and the standard
% error of the mean
Phi = [Phireal,mean(ComCount(:,6)), std(ComCount(:,6))/sqrt(numberOfFrames)];


%6. Store data for plotting
%{
FileExists = exist('Polang.mat');
if FileExists == 2
    save ('Polang.mat',Phi45,'-append');
else
    save ('Polang.mat',['Phi',int2str(Phireal)]);
end
%}
if flagInitialRun
        PhiPlot = Phi;
else
        PhiPlot = cat(1, PhiPlot, Phi);
end

