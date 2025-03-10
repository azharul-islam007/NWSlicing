function [] = setup_network_model_dl_f()

commSupportPackageCheck("CST_WINNER2");

s = rng(21); % For repeatability

az = -180:179; % 1-degree spacing
pattern = cat(1,shiftdim(winner2.dipole(az,45),-1), shiftdim(winner2.dipole(az,-45),-1), shiftdim(winner2.dipole(az,45),-1), shiftdim(winner2.dipole(az,-45),-1));
ULA4_Bs = winner2.AntennaArray('ULA',4,0.5, 'FP-ECS',pattern,'Azimuth',az);

pattern = cat(1,shiftdim(winner2.dipole(az,45),-1), shiftdim(winner2.dipole(az,-45),-1));
ULA2_Ue_Car = winner2.AntennaArray('ULA',2,0.5, 'FP-ECS',pattern,'Azimuth',az);

AA(1) = ULA4_Bs;
AA(2) = ULA2_Ue_Car;

M_Bs = 1;
M_Ue = 4;
M_car_cellular_dl = 2;
M_car_sidelink = 2;
M_car_others = 2;

TxAntArrayIdx = {1 * ones(1, M_Bs); 2 * ones(1, M_car_sidelink)}; % 2 Index in antenna array inventory vector
RxAntArrayIdx = 2 * ones(1, M_Ue + M_car_cellular_dl + M_car_others);     % 6 Index in antenna array inventory vector
numLinks = M_Ue + M_car_cellular_dl + M_car_sidelink * M_car_others;            % 8 Number of links
radius = 500;
range = 2 * radius;             % Layout range (meters)
cfgLayout = winner2.layoutparset(RxAntArrayIdx, TxAntArrayIdx, numLinks, AA, range);

numTx = M_Bs + M_car_sidelink;
numRx = length(RxAntArrayIdx);
numLink = M_Ue + M_car_cellular_dl + M_car_sidelink * M_car_others;

firstTx_Bs = 1;
firstTx_car_sidelink = firstTx_Bs + 1;
lastTx_car_sidelink = firstTx_Bs + M_car_sidelink;

firstRx_Ue = numTx + 1;
lastRx_Ue = numTx + M_Ue;
firstRx_car_cellular = lastRx_Ue + 1;
lastRx_car_cellular = lastRx_Ue + M_car_cellular_dl;
firstRx_car_others = lastRx_car_cellular + 1;
lastRx_car_other = lastRx_car_cellular + M_car_others;

pairing_Bs_Ue = [ones(1, M_Ue); firstRx_Ue : lastRx_Ue];
pairing_Bs_car_cellular = [ones(1, M_car_cellular_dl); firstRx_car_cellular : lastRx_car_cellular];
pairing_sidelink = [repelem(firstTx_car_sidelink: lastTx_car_sidelink, M_car_others); 
    repmat(firstRx_car_others: lastRx_car_other, 1, M_car_others)];

cfgLayout.Pairing = [pairing_Bs_Ue pairing_Bs_car_cellular pairing_sidelink];

% cfgLayout.Pairing = [1 1 1 1 1 1 2 2; 
%                      3 4 5 6 7 8 9 10];  % Index in cfgLayout.Stations
cfgLayout.ScenarioVector = 3 * ones(1, numLink); % [3 3 3 3 3 3 3 3];     % 3 for B1
cfgLayout.PropagConditionVector = [1 1 1 1 1 1 1 1 1 1];  % 0 for NLOS

% Set up positions for BS sectors. 
sqrt3 = sqrt(3);

cfgLayout.Stations(1).Pos(1:2) = [500;  round(radius / 2 * sqrt3)];
% Set up positions for cars transmitting in sidelink mode 
first_station_car_sidelink = M_Bs + 1;
last_station_car_sidelink = M_Bs + M_car_sidelink;
for index = first_station_car_sidelink : last_station_car_sidelink
    offset = index - first_station_car_sidelink + 1;
    cfgLayout.Stations(index).Pos(1:2) = [500 - offset * 100;  round(radius / 2 * sqrt3 + 30 * offset)]; % 50m from 1st BS
end

% Set up MS positions
first_station_UE = M_Bs + M_car_sidelink + 1;
last_station_UE = M_Bs + M_car_sidelink + M_Ue;
for index = first_station_UE : last_station_UE
    offset = index - first_station_UE + 1;
    cfgLayout.Stations(index).Pos(1:2) = [50 * offset; 50 * offset]; % 55m and 65m from 2nd and 3rd BSs respectively
end

first_station_car_cellular = M_Bs + M_car_sidelink + M_Ue + 1;
last_station_car_cellular = M_Bs + M_car_sidelink + M_Ue + M_car_cellular_dl;
for index = first_station_car_cellular : last_station_car_cellular
    offset = index - first_station_car_cellular + 1;
    cfgLayout.Stations(index).Pos(1:2) = [500 + offset  * 100;  round(radius / 2 * sqrt3 - 10 * offset)]; % 55m and 65m from 2nd and 3rd BSs respectively
end

first_station_car_others = M_Bs + M_car_sidelink + M_Ue + M_car_cellular_dl + 1;
last_station_car_others =  M_Bs + M_car_sidelink + M_Ue + M_car_cellular_dl + M_car_others;
for index = first_station_car_others : last_station_car_others
    offset = index - first_station_car_others + 1;
    cfgLayout.Stations(index).Pos(1:2) = [offset * 100;  round(radius / 2 * sqrt3 - 10 * offset)]; % 55m and 65m from 2nd and 3rd BSs respectively
end

% Randomly draw MS velocity
for i = numTx + (1:numRx)
    cfgLayout.Stations(i).Velocity = rand(3,1) - 0.5;
end

% Get all BS sector and MS positions
BSPos = cell2mat({cfgLayout.Stations(1:numTx).Pos});
MSPos = cell2mat({cfgLayout.Stations(numTx+1:end).Pos});

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
xlim([0 1000]); ylim([0 1000]);
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