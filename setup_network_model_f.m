function [] = setup_network_model_f()

commSupportPackageCheck("CST_WINNER2");

s = rng(21); % For repeatability

az = -180:179; % 1-degree spacing
pattern = cat(1,shiftdim(winner2.dipole(az,45),-1), shiftdim(winner2.dipole(az,-45),-1), shiftdim(winner2.dipole(az,45),-1), shiftdim(winner2.dipole(az,-45),-1));
ULA4_Bs = winner2.AntennaArray('ULA',4,0.5, 'FP-ECS',pattern,'Azimuth',az);

pattern = cat(1,shiftdim(winner2.dipole(az,45),-1), shiftdim(winner2.dipole(az,-45),-1));
ULA2_Ue_Car = winner2.AntennaArray('ULA',2,0.5, 'FP-ECS',pattern,'Azimuth',az);

AA(1) = ULA4_Bs;
AA(2) = ULA2_Ue_Car;

M_Ue = 4;
M_car_cellular_dl = 2;
M_car_sidelink = 1;
M_car_others = 2;

BSIdx = {1; 2 * ones(M_car_sidelink, 1)}; % 2 Index in antenna array inventory vector
MSIdx = 2 * ones(1, M_Ue + M_car_cellular_dl + M_car_others);     % 6 Index in antenna array inventory vector
numLinks = M_Ue + M_car_cellular_dl + M_car_sidelink * M_car_others;            % 8 Number of links
range = 500;             % Layout range (meters)
cfgLayout = winner2.layoutparset(MSIdx, BSIdx, numLinks, AA, range);

cfgLayout.Pairing = [1 1 1 1 1 1 2 2; 
                     3 4 5 6 7 8 9 10];  % Index in cfgLayout.Stations
cfgLayout.ScenarioVector = [3 3 3 3 3 3 3 3];     % 3 for B1
cfgLayout.PropagConditionVector = [0 0 0 0 0 0 0 0];  % 0 for NLOS

% Number of BS sectors and MSs in the system
numBSSect = sum(cfgLayout.NofSect);
numMS = length(MSIdx);

% Set up positions for BS sectors. Same position for the 
% third, fourth and fifth sectors as they belong to one BS.
cfgLayout.Stations(1).Pos(1:2) = [50;  150];

% Set up MS positions
cfgLayout.Stations(2).Pos(1:2) = [10;  180]; % 50m from 1st BS
cfgLayout.Stations(3).Pos(1:2) = [60;  50];  % 111.8m from 1st BS
cfgLayout.Stations(4).Pos(1:2) = [194; 117]; % 55m and 65m from 2nd and 3rd BSs respectively
cfgLayout.Stations(5).Pos(1:2) = [260; 270]; % 120.4m from 3rd BS
cfgLayout.Stations(6).Pos(1:2) = [295; 90]; % 75m from 3rd BS
cfgLayout.Stations(7).Pos(1:2) = [295; 90]; % 75m from 3rd BS
cfgLayout.Stations(8).Pos(1:2) = [295; 110]; % 75m from 3rd BS
cfgLayout.Stations(9).Pos(1:2) = [295; 130]; % 75m from 3rd BS
cfgLayout.Stations(10).Pos(1:2) = [295; 150]; % 75m from 3rd BS

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
cfgWim.CenterFrequency = 2.6e9;
cfgWim.UniformTimeSampling = "no";
cfgWim.ShadowingModelUsed = "yes";
cfgWim.PathLossModelUsed = "yes";
cfgWim.RandomSeed = 31415926;       % For repeatability

WINNERChan = comm.WINNER2Channel(cfgWim,cfgLayout);
chanInfo = info(WINNERChan);

txSig = cellfun(@(x) [ones(1,x);zeros(frameLen-1,x)], ...
    num2cell(chanInfo.NumBSElements)',UniformOutput=false);

[rxSig, path_gains] = WINNERChan(txSig); % Pass impulse signal through each link

path_gains









end