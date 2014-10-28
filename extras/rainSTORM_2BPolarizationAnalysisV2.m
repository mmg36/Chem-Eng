%rainSTORM_2BPolarizationAnalysis
%This script uses the main RainSTORM Software to analyse 2 stacks of x-y
%polarized images

% Made by Daan van Kleef & Mehdi Goudarzi


% 1. Use rainSTORM to analyse raw image data
%    Then read the data from the Matlab struct into a useful array...
clear
startup  % Run rainSTORM once (manually), on the X-polarisation data
uiwait 

% Read in as an array!
Nxtable = struct2table((params.localization.results.SupResParams));
for count = 1:size(Nxtable,2)
    Nxdummy = double(table2array(Nxtable(:,count)));
    if count == 1
        Nxarray = Nxdummy
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
        Nyarray = Nydummy
    else
        Nyarray = cat(2,Nyarray,Nydummy);
    end
end


%Nytable = double(table2array(Nytable));
%Ny =double(table2array(Nytable(:,5)));

% 2. Identify matching images, and force search for paired signal if this
%    is not yet found
%    Method 1: assume fluorophore position is known (from simulation
%    params)

% Method 1:
NxtableCorrect = Nxtable;
NxDistance     = sqrt( (Nxtable(:,2)-32).^2 + (Nxtable(:,3)-30).^2 );
NxtableCorrect(NxDistance > 1, :) = [];

NytableCorrect = Nytable;
NyDistance     = sqrt( (Nytable(:,2)-32).^2 + (Nytable(:,3)-30).^2 );
NytableCorrect(NyDistance > 1, :) = [];

numberOfFrames = params.localization.results.numberOfFrames;

NxTC_ytableMatch = -ones(size(NxtableCorrect,1),1);

% 3. Evaluate polarisation (or orientation) of observed molecules:
