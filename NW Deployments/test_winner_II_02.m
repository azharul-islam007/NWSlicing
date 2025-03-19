commSupportPackageCheck("CST_WINNER2");

s = rng(21); % For repeatability

AA(1) = winner2.AntennaArray("UCA",16,0.3);
AA(2) = winner2.AntennaArray("UCA",12,0.3);
AA(3) = winner2.AntennaArray("UCA",8,0.3);
AA(4) = winner2.AntennaArray("UCA",4,0.05);

BSIdx = {1; 2; [3 3 3]}; % Index in antenna array inventory vector
MSIdx = [4 4 4 4 4];     % Index in antenna array inventory vector
numLinks = 7;            % Number of links
range = 300;             % Layout range (meters)
cfgLayout = winner2.layoutparset(MSIdx,BSIdx,numLinks,AA,range);

cfgLayout.Pairing = [1 1 2 3 3 4 6; 6 7 8 8 9 10 7];  % Index in cfgLayout.Stations
cfgLayout.ScenarioVector = [6 6 13 13 11 11 5];     % 6 for B4, 11 for C2 and 13 for C4
cfgLayout.PropagConditionVector = [0 0 0 0 0 0 0];  % 0 for NLOS

% Number of BS sectors and MSs in the system
numBSSect = sum(cfgLayout.NofSect);
numMS = length(MSIdx);

% Set up positions for BS sectors. Same position for the 
% third, fourth and fifth sectors as they belong to one BS.
cfgLayout.Stations(1).Pos(1:2) = [50;  150];
cfgLayout.Stations(2).Pos(1:2) = [150; 150];
cfgLayout.Stations(3).Pos(1:2) = [250; 150];
cfgLayout.Stations(4).Pos(1:2) = [250; 150];
cfgLayout.Stations(5).Pos(1:2) = [250; 150];

% Set up MS positions
cfgLayout.Stations(6).Pos(1:2) = [10;  180]; % 50m from 1st BS
cfgLayout.Stations(7).Pos(1:2) = [60;  50];  % 111.8m from 1st BS
cfgLayout.Stations(8).Pos(1:2) = [194; 117]; % 55m and 65m from 2nd and 3rd BSs respectively
cfgLayout.Stations(9).Pos(1:2) = [260; 270]; % 120.4m from 3rd BS
cfgLayout.Stations(10).Pos(1:2) = [295; 90]; % 75m from 3rd BS

% Randomly draw MS velocity
for i = numBSSect + (1:numMS)
    cfgLayout.Stations(i).Velocity = rand(3,1) - 0.5;
end

% Get all BS sector and MS positions
BSPos = cell2mat({cfgLayout.Stations(1:numBSSect).Pos});
MSPos = cell2mat({cfgLayout.Stations(numBSSect+1:end).Pos});

scrsz = get(groot,"ScreenSize");
figSize = min(scrsz([3,4]))/2.3;
figure( ...
    Position=[scrsz(3)*.5-figSize/2,scrsz(4)*.7-figSize/2,figSize,figSize]);
hold on; 
grid on;
hBS = plot(BSPos(1,:),BSPos(2,:),"or");   % Plot BS
hMS = plot(MSPos(1,:),MSPos(2,:),"xb");   % Plot MS
for linkIdx = 1:numLinks                  % Plot links
    pairStn = cfgLayout.Pairing(:,linkIdx);
    pairPos = cell2mat({cfgLayout.Stations(pairStn).Pos});
    plot(pairPos(1,:),pairPos(2,:),"-b");
end
xlim([0 300]); ylim([0 300]);
xlabel("X Position (meters)");
ylabel("Y Position (meters)")
legend([hBS, hMS],"BS","MS",location="northwest");

frameLen = 1600; % Number of samples to be generated

cfgWim = winner2.wimparset;
cfgWim.NumTimeSamples = frameLen;
cfgWim.IntraClusterDsUsed = "yes";
cfgWim.CenterFrequency = 5.25e9;
cfgWim.UniformTimeSampling = "no";
cfgWim.ShadowingModelUsed = "yes";
cfgWim.PathLossModelUsed = "yes";
cfgWim.RandomSeed = 31415926;       % For repeatability

