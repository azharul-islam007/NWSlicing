classdef Network_Model_c < handle
    properties
        gNB_Bs; % a unique gNB base station
        eMBB_UEs = {};  % a list of eMBB UEs in the network
        vehicular_UEs = {}; % a list of vehicular UEs in the network
        mMTC_UEs = {} % a list of mMTC UEs in the network
        M_eMBB; % number of eMBB UEs
        M_veh; % number of vehicular UEs
        M_mMTC; % number of mMTC UEs
        bsId; % Base station Id, fixed to be 1
        radius_min; % eMBB UEs or mMTC UEs are randomly distrubted at a minimum distance away from the gNB base station
        radius_max; % the maximum cell coverage radius away from the gNB base station
        
        pd_radius; % a uniform distribution for the UEs positions away from the gNB
        pd_theta; % a uniform distribution for the UEs around the center of the gNB
        pd_theta1;
        pd_veh_ue_direction; % vehicular UEs are randomly driving from left to right or right to left
        pd_veh_ue_x_random; % the uniform distribution of the initial positions of vehicular UEs at the start of the simulation
        pd_veh_ue_x_l2r; % the uniform distribution of the position (x) for vehicular UEs when they arrive from left to right
        pd_veh_ue_x_r2l; % the uniform distribution of the position (x) for vehicular UEs when they arrive from right to left
        direction_str; % = {'l2r', 'r2l'}

        next_ueId; % keep track of the next ueId, used when new UEs arrive
        pd_eMBB_session_N; % the Poission distribution of the generation of the eMBB sessions
        pd_eMBB_session_T; % the uniform distribution of the duration of an eMBB session
        pd_eMBB_UE_birth_N; % the arriving distriubtion of the eMBB UE

        pd_veh_packet_N; % the Poission distribution of the generation the V2X packets
        pd_veh_UE_birth_N; % the arriving distribution of the vehicular UEs (cars)

        lambda_e; % arriving rate of the eMBB session 
        lambda_m; % arrivate rate of the eMBB UE

        lambda_niu; % packet arriving rate of V2X service
        lambda_a; % vehicular UEs arriving rate

        sim_time; % number of simulation drops

        winner2_channel_smp_time; % number of samples to be generated when generate WINNER2 impulse channel
        WINNERChan; % the generated WINNER2 channel, which should be updated in every simulation drop
        ULA4_Bs; % ULA antenna array with 4 elements, used in gNB, with antenna gain approximately 5dB
        ULA2_Ue_Car; % ULA antenna array with 2 element, used in all the UEs (eMBB, mMTC, or Vehicle UEs)
    end

    methods
        function self = Network_Model_c(M_eMBB_i, M_veh_i, M_mMTC_i, lambda_e_i, sim_time_i)
            % Construct and initialize the network models and its
            % properties (as descripted in the class properties)
            % Refer to Table 1 (simulation parameters) in the paper
            s = rng(21); % For repeatability
            self.sim_time = sim_time_i;
            self.M_eMBB = M_eMBB_i;
            self.M_veh = M_veh_i;
            self.M_mMTC = M_mMTC_i;
            self.lambda_e = lambda_e_i; % varied for performance comparison
            self.lambda_m = 1; % 1 UE/second
            self.lambda_niu = 1; % 1 packet/second
            self.lambda_a = 1; % 1 UE/second
            self.next_ueId = 1; 
            self.radius_min = 25;
            self.radius_max = Constants_c.radius / 2 * sqrt(3);
            self.pd_radius = makedist("Uniform", "lower", self.radius_min, "upper", self.radius_max); % uniform distribution
            self.pd_theta = makedist("Uniform", "lower", pi/18, "upper", 17*pi/18); % uniform distribution
            self.pd_theta1 = makedist("Uniform", "lower", 19*pi/18, "upper", 35*pi/18); % uniform distribution
            self.pd_veh_ue_direction = makedist("Binomial", "N", 1, "p", 0.5); % 0 - left2right; 1 - right2left
            self.pd_veh_ue_x_random = makedist("Uniform", "lower", 0, "upper", Constants_c.diameter); % uniform distribution
            self.direction_str = {'l2r', 'r2l'};
            v = Constants_c.velocity;
            t = Constants_c.T_drop;
            delta_d = v * t; % a distance for a vehicle moving forward within on simulation drop
            self.pd_veh_ue_x_l2r = makedist("Uniform", "lower", 0, "upper", delta_d); % uniform distribution
            self.pd_veh_ue_x_r2l = makedist("Uniform", ... % uniform distribution
                "lower", Constants_c.diameter, ...
                "upper", Constants_c.diameter + delta_d);

            lambda_e_per_drop = self.lambda_e / (1 / Constants_c.T_drop); 
            self.pd_eMBB_session_N = makedist("Poisson", "lambda", lambda_e_per_drop); % Poission distribution
            self.pd_eMBB_session_T = makedist("Uniform", "lower", Constants_c.T_session_eMBB_min, "upper", Constants_c.T_session_eMBB_max); % uniform distribution
            lambda_m_per_drop = self.lambda_m / (1 / Constants_c.T_drop);
            self.pd_eMBB_UE_birth_N = makedist("Poisson", "lambda", lambda_m_per_drop); % Poission distribution

            lambda_niu_per_drop = self.lambda_niu / (1 / Constants_c.T_drop);
            self.pd_veh_packet_N = makedist("Poisson", "lambda", lambda_niu_per_drop); % Poission distribution
            lambda_a_per_drop = self.lambda_a / (1 / Constants_c.T_drop);
            self.pd_veh_UE_birth_N = makedist("Poisson", "lambda", lambda_a_per_drop);    % Poission distribution

            self.winner2_channel_smp_time = 1600; % recommeded value from mathworld's default configuration
        end

        function [] = init_winnerII_model(self)
            % initialize the WINNER2 model according to this example:
            % https://ww2.mathworks.cn/help/comm/ref/winner2.antennaarray.html
            % https://ww2.mathworks.cn/help/comm/ug/simultaneous-simulation-of-multiple-fading-channels-with-winner-ii-channel-model.html
            commSupportPackageCheck("CST_WINNER2");      
            s = rng(21); % For repeatability
            az = -180:179; % 1-degree spacing
            pattern = cat(1, shiftdim(winner2.dipole(az, 45), -1), ...
                shiftdim(winner2.dipole(az, -45), -1), ...
                shiftdim(winner2.dipole(az, 45), -1), ...
                shiftdim(winner2.dipole(az, -45), -1));
            self.ULA4_Bs = winner2.AntennaArray('ULA', 4, 0.5, 'FP-ECS', pattern, 'Azimuth', az);
            
            pattern = cat(1,shiftdim(winner2.dipole(az, 45), -1), ...
                shiftdim(winner2.dipole(az, -45), -1));
            self.ULA2_Ue_Car = winner2.AntennaArray('ULA', 2, 0.5, 'FP-ECS', pattern, 'Azimuth',az);
        end

        function [] = init_simulation(self)
            % initialize the positions of the gNB, eMBB, mMTC, vehicular
            % UEs
            self.bsId = 1; % fixed bsId = 1
            self.gNB_Bs = BS_c(self.bsId, Constants_c.gNB_Bs_pos_x, ...
                Constants_c.gNB_Bs_pos_y, Constants_c.gNB_Bs_height);
            self.eMBB_UEs = cell(self.M_eMBB, 1);
            % eMBB UEs are uniformly distribution around the base station
            % with a minimum and maximum distrance
            for count = 1 : self.M_eMBB
                r = self.pd_radius.random();
                theta = self.pd_theta.random();
                x = self.gNB_Bs.pos_x + r * cos(theta);
                y = 481 + r * sin(theta);
                self.eMBB_UEs{count} = eMBB_UE_c(self.next_ueId, x, y, ...
                    Constants_c.UE_height);
                self.next_ueId = self.next_ueId + 1;
            end
            % mMTC UEs are uniformly distribution around the base station
            % with a minimum and maximum distrance
            self.mMTC_UEs = cell(self.M_mMTC, 1);
            for count = 1 : self.M_mMTC
                r = self.pd_radius.random();
                theta1 = self.pd_theta1.random();
                x = self.gNB_Bs.pos_x + r * cos(theta1);
                y = 385 + r * sin(theta1);
                T_period = 20; % assuming the eMTCs UEs are sending packets every T_period
                T_offset = randsample(1: T_period, 1);
                packet_size = 100; % assuming the eMTC UEs are sending a packet of size 100 bytes
                p_Tx = 0.5; % a probability that the eMTC UEs will wake up to send a packet
                self.mMTC_UEs{count} = mMTC_UE_c(self.next_ueId, x, y, ...
                    Constants_c.UE_height, T_period, T_offset, packet_size, p_Tx);
                self.next_ueId = self.next_ueId + 1;                
            end
            % vehicular UEs are randomly distrubtion in the high ways,
            % either driving from left to right or from right to left
            self.vehicular_UEs = cell(self.M_veh, 1);
            for count = 1 : self.M_veh
                rn = self.pd_veh_ue_direction.random(); % rn = 0 or 1
                direction = self.direction_str{rn + 1};
                if strcmpi(direction, 'l2r') % from left to right
                    x = self.pd_veh_ue_x_random.random();
                    lane = randsample(1: Constants_c.nr_lane_one_direction, 1);
                    y = Constants_c.gNB_Bs_pos_y - Constants_c.d_gNB_highway - (lane - 0.5) * Constants_c.d_lane;
                    self.vehicular_UEs{count} = vehicular_UE_c(self.next_ueId, x, y, ...
                        Constants_c.UE_height, "l2r");
                    self.next_ueId = self.next_ueId + 1;
                elseif strcmpi(direction, 'r2l')
                    x = self.pd_veh_ue_x_random.random();
                    lane = randsample(Constants_c.nr_lane_one_direction + 1: Constants_c.nr_lane_bi_direction, 1);
                    y = Constants_c.gNB_Bs_pos_y - Constants_c.d_gNB_highway - (lane - 0.5) * Constants_c.d_lane;
                    self.vehicular_UEs{count} = vehicular_UE_c(self.next_ueId, x, y, ...
                        Constants_c.UE_height, "r2l");
                    self.next_ueId = self.next_ueId + 1;                    
                else
                    error('either 12r or r2l is allowed....')
                end
            end
        end

        function [] = prepare_simulation_drop(self, t_i)
            % before the start of every simulation drop, there may be new
            % eMBB UEs arrived and replace some old eMBB UEs; and there may
            % be vehicular UEs arriving either driving from left to right
            % or from right to left. Those vehicular UEs driving out of the
            % highway will be removed.
            % For every eMBB UEs, there may by new sessions generated or
            % old sessions died.
            % For every vehicular UEs, there may be new packets arrived for
            % transmission immediately.

            % create a number of new eMBB UE (L_eMBB) with the probability of arrival
            % rate λm
            L_eMBB = min(self.pd_eMBB_UE_birth_N.random(), self.M_eMBB);
            if L_eMBB > 0
                % Replace randomly L_eMBB UEs in the network at the same time
                eMBB_UE_indices = randsample(self.M_eMBB, L_eMBB);
                for index = eMBB_UE_indices
                    r = self.pd_radius.random();
                    theta = self.pd_theta.random();
                    x = self.gNB_Bs.pos_x + r * cos(theta);
                    y = self.gNB_Bs.pos_y + r * sin(theta);
                    self.eMBB_UEs{index} = eMBB_UE_c(...
                        self.next_ueId, x, y, ...
                        Constants_c.UE_height);
                    self.next_ueId = self.next_ueId + 1;
                end
            end

            % For each eMBB UE, it should first refresh the session info
            % before the start of each simulation drop
            for m = 1 : self.M_eMBB
                self.eMBB_UEs{m}.refresh_session(t_i); % first
            end

            % For each eMBB UE, it should trigger a number of K_(m,eMBB)  new sessions
            for m = 1 : self.M_eMBB
                K_m_eMBB = self.pd_eMBB_session_N.random();
                if K_m_eMBB > 0
                    for count = 1 : K_m_eMBB
                        T_m_eMBB = self.pd_eMBB_session_T.random();
                        self.eMBB_UEs{m}.add_session(t_i, T_m_eMBB);
                    end
                end
            end
            
            % Remove the vehicular UE that has driven out of the coverage
            for count = 1 : self.M_veh
                self.vehicular_UEs{count}.move_forward_by_one_drop();
                if self.vehicular_UEs{count}.is_out_of_highway()
                    self.vehicular_UEs{count} = [];
                end
            end
            self.vehicular_UEs = self.vehicular_UEs(~cellfun('isempty', self.vehicular_UEs));
            self.M_veh = numel(self.vehicular_UEs);

            % Generate a number of new vehicular UE 𝐿𝑉2 with the probability 𝑃(𝐿𝑉2𝑋 = 𝑙) distribution 
            % conforming to the Poisson Distribution 𝑃(𝜆𝑎)
            L_V2X = self.pd_veh_UE_birth_N.random();
            if L_V2X > 0
                for count = 1 : L_V2X
                    rn = self.pd_veh_ue_direction.random(); % rn = 0 or 1
                    direction = self.direction_str{rn + 1};
                    if strcmpi(direction, 'l2r') % from left to right
                        x = self.pd_veh_ue_x_l2r.random();
                        lane = randsample(1: Constants_c.nr_lane_one_direction, 1);
                        y = Constants_c.gNB_Bs_pos_y - Constants_c.d_gNB_highway - (lane - 0.5) * Constants_c.d_lane;
                        self.vehicular_UEs{self.M_veh + count} = vehicular_UE_c(self.next_ueId, x, y, ...
                            Constants_c.UE_height, "l2r");
                        self.next_ueId = self.next_ueId + 1;
                    elseif strcmpi(direction, 'r2l')
                        x = self.pd_veh_ue_x_r2l.random();
                        lane = randsample(Constants_c.nr_lane_one_direction + 1: Constants_c.nr_lane_bi_direction, 1);
                        y = Constants_c.gNB_Bs_pos_y - Constants_c.d_gNB_highway - (lane - 0.5) * Constants_c.d_lane;
                        self.vehicular_UEs{self.M_veh + count} = vehicular_UE_c(self.next_ueId, x, y, ...
                            Constants_c.UE_height, "r2l");
                        self.next_ueId = self.next_ueId + 1;                    
                    else
                        error('either 12r or r2l is allowed....')
                    end
                end
                self.M_veh = self.M_veh + L_V2X;
            end
            
            % For each vehicular UE, it should trigger a number of 𝐾𝑚,𝑉2𝑋 new packets with the probability 
            % 𝑃(𝐾𝑚,𝑉2𝑋 = 𝑘) distribution conforming to the Poisson Distribution 𝑃(𝜆𝑣)
            for count = 1 : self.M_veh
                K_m_V2X = self.pd_veh_packet_N.random();
                if K_m_V2X > 0
                    self.vehicular_UEs{count}.add_packets(K_m_V2X);
                end
            end
        end       


        function [h_gNB_Bs_o, h_eMBB_UE_o, h_mMTC_UE_o] = plot(self)
            % plot the network region,  the highway and the UEs
            scrsz = get(groot,"ScreenSize");
            figSize = min(scrsz([3,4]))/2.3;
            figure(...
                Position=[scrsz(3)*.5-figSize/2,scrsz(4)*.7-figSize/2,figSize,figSize]);            
            hold on; 
            grid on;
            % plot the hexagonal region and the highway lanes
            r = Constants_c.radius;
            r_0_5 = r * 0.5;
            r_0_86 = r * sqrt(3) * 0.5;
            x1 = [0, r_0_86];
            x2 = [r_0_5, 0];
            x3 = [r + r_0_5, 0];
            x4 = [2 * r, r_0_86];
            x5 = [r + r_0_5, 2 * r_0_86];
            x6 = [r_0_5, 2 * r_0_86];
            plot([x1(1), x2(1)], [x1(2), x2(2)], '-k');
            plot([x2(1), x3(1)], [x2(2), x3(2)], '-k');
            plot([x3(1), x4(1)], [x3(2), x4(2)], '-k');
            plot([x4(1), x5(1)], [x4(2), x5(2)], '-k');
            plot([x5(1), x6(1)], [x5(2), x6(2)], '-k');
            plot([x6(1), x1(1)], [x6(2), x1(2)], '-k');
            lane_offset = Constants_c.gNB_Bs_pos_y - Constants_c.d_gNB_highway;
            w = Constants_c.d_lane;
            for lane = 0 : Constants_c.nr_lane_bi_direction
                left = [0, lane_offset - (lane * w)]; right = [Constants_c.diameter, lane_offset - (lane * w)];
                if lane==3
                    plot([left(1), right(1)], [left(2), right(2)], '-r');
                else
                    plot([left(1), right(1)], [left(2), right(2)], '-k');
                end
            end
            annotation('textarrow',[0.2 0.25],[0.47 0.47],'LineStyle','-','color',[1 0 1]);
            annotation('textarrow',[0.2 0.25],[0.482 0.482],'LineStyle','-','color',[1 0 1]);
            annotation('textarrow',[0.2 0.25],[0.495 0.495],'LineStyle','-','color',[1 0 1]);
            annotation('textarrow',[0.25 0.2],[0.457 0.457],'LineStyle','-','color',[0 1 0]);
            annotation('textarrow',[0.25 0.2],[0.444 0.444],'LineStyle','-','color',[0 1 0]);
            annotation('textarrow',[0.25 0.2],[0.431 0.431],'LineStyle','-','color',[0 1 0]);
          
            % plot gNB_BS
            h_gNB_Bs_o = plot(self.gNB_Bs.pos_x, self.gNB_Bs.pos_y, "or");
            % plot eMBB UE
            eMBB_UE_pos_x = zeros(self.M_eMBB, 1);
            eMBB_UE_pos_y = zeros(self.M_eMBB, 1);
            for count = 1 : self.M_eMBB
                eMBB_UE_pos_x(count) = self.eMBB_UEs{count}.pos_x;
                eMBB_UE_pos_y(count) = self.eMBB_UEs{count}.pos_y;
            end
            h_eMBB_UE_o = plot(eMBB_UE_pos_x, eMBB_UE_pos_y, "*b");

            % plot mMTC UE
            mMTC_UE_pos_x = zeros(self.M_mMTC, 1);
            mMTC_UE_pos_y = zeros(self.M_mMTC, 1);
            for count = 1 : self.M_mMTC
                mMTC_UE_pos_x(count) = self.mMTC_UEs{count}.pos_x;
                mMTC_UE_pos_y(count) = self.mMTC_UEs{count}.pos_y;
            end
            h_mMTC_UE_o = plot(mMTC_UE_pos_x, mMTC_UE_pos_y, "LineStyle", "none", "Marker", "square", "Color", "#A2142F");            

            % plot the links between gNB and eMBB UEs/mMTC UEs
            for count = 1 : self.M_eMBB
                x = [self.gNB_Bs.pos_x self.eMBB_UEs{count}.pos_x];
                y = [self.gNB_Bs.pos_y self.eMBB_UEs{count}.pos_y];
                plot(x, y, "--b");
                text(x(2),y(2), sprintf('sinr=%.1f dB', self.eMBB_UEs{count}.sinr_dB));
            end
            for count = 1 : self.M_mMTC
                x = [self.gNB_Bs.pos_x self.mMTC_UEs{count}.pos_x];
                y = [self.gNB_Bs.pos_y self.mMTC_UEs{count}.pos_y];
                plot(x, y, "--g");
                text(x(2),y(2), sprintf('sinr=%.1f dB', self.mMTC_UEs{count}.sinr_dB));
            end            
        end

        function [] = run_winner2_model_simulation_drop(self, t_i)
            % should be override by its subclass
        end

        function [] = start_simulation(self)
            self.init_simulation();
            self.init_winnerII_model();
            for t = 1 : self.sim_time         % sim_time: num. of simulation drop=10
                fprintf(1, '.');
                self.prepare_simulation_drop(t);                
                self.run_winner2_model_simulation_drop(t);
                if (mod(t-1, floor(1/Constants_c.T_drop)) == 0) || (t == self.sim_time) % when t=1 or t=10, plot starts
                    self.plot();
                end                
            end
        end
    end

end