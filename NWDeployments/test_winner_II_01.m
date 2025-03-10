BSAA  = winner2.AntennaArray('UCA', 8, 0.02);  % UCA-8 array for BS
MSAA1 = winner2.AntennaArray('ULA', 2, 0.01);  % ULA-2 array for MS
MSAA2 = winner2.AntennaArray('ULA', 4, 0.005); % ULA-4 array for MS

MSIdx = [2 3]; 
BSIdx = {1}; 
K = 2; 
rndSeed = 5;
cfgLayout = winner2.layoutparset(MSIdx,BSIdx, ...
    K,[BSAA,MSAA1,MSAA2],[],rndSeed);

BSPos  = cfgLayout.Stations(cfgLayout.Pairing(1,1)).Pos;
MS1Pos = cfgLayout.Stations(cfgLayout.Pairing(2,1)).Pos;
MS2Pos = cfgLayout.Stations(cfgLayout.Pairing(2,2)).Pos;

plot3(BSPos(1),BSPos(2),BSPos(3),'bo', ...
    MS1Pos(1),MS1Pos(2),MS1Pos(3),'rs', ...
    MS2Pos(1),MS2Pos(2),MS2Pos(3),'rd');
grid on;
xlim([0 500]);
ylim([0 500]);
zlim([0 35]);
xlabel('X-position (m)');
ylabel('Y-position (m)');
zlabel('Elevation (m)');
legend('BS','MS1','MS2','Location','northeast');