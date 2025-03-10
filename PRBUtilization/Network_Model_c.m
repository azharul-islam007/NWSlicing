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

        % the distribution of the index i \in {1, 2, ...,
        % log2(mMTC_max_Packet_Size)}, the actual pakcet size is 2^i bytes
        pd_mMTC_log2_packet_size; 
        % the dsitribution of the transmission probability (i.e. the
        % probability of an mMTC UE waking up to perform transmission)
        pd_mMTC_tx;

        sim_nr_drops_per_round; % number of simulation drops per round

        winner2_channel_smp_time; % number of samples to be generated when generate WINNER2 impulse channel
        WINNERChan; % the generated WINNER2 channel, which should be updated in every simulation drop
        ULA4_Bs; % ULA antenna array with 4 elements, used in gNB, with antenna gain approximately 5dB
        ULA2_Ue_Car; % ULA antenna array with 2 element, used in all the UEs (eMBB, mMTS, or Vehicle UEs)

        nw_characteristics_v2x;
        nw_characteristics_eMBB;
        nw_characteristics_mMTC;

        alpha_v2x;  % slicing ratio for V2X
        alpha_eMBB; % slicing ratio for eMBB
        alpha_mMTC; % slicing ratio for mMTC

        prb_utilization_ratio_v2x;
        prb_utilization_ratio_eMBB;
        prb_utilization_ratio_mMTC;
        prb_utilization_ratio_total;
    end

    methods
        function self = Network_Model_c(M_eMBB_i, M_veh_i, M_mMTC_i, lambda_e_i, sim_nr_drops_per_round_i)
            % Construct and initialize the network models and its
            % properties (as descripted in the class properties)
            % Refer to Table 1 (simulation parameters) in the paper
            s = rng(21); % For repeatability
            self.sim_nr_drops_per_round = sim_nr_drops_per_round_i;
            self.M_eMBB = M_eMBB_i;
            self.M_veh = M_veh_i;
            self.M_mMTC = M_mMTC_i;
            self.lambda_e = lambda_e_i; % varified for performance comparison
            self.lambda_m = Constants_c.lambda_m; 
            self.lambda_niu = Constants_c.lambda_niu; 
            self.lambda_a = Constants_c.lambda_a;            
            self.radius_min = Constants_c.radius_min;
            self.radius_max = Constants_c.radius_max;

            self.pd_radius = makedist("Uniform", "lower", self.radius_min, "upper", self.radius_max); % uniform distribution
            self.pd_theta = makedist("Uniform", "lower", 0, "upper", 2 * pi); % uniform distribution
            self.pd_veh_ue_direction = makedist("Binomial", "N", 1, "p", 0.5); % 0 - left2right; 1 - right2left; vehicles are driving left or right equally
            self.pd_veh_ue_x_random = makedist("Uniform", "lower", 0, "upper", Constants_c.diameter); % uniform distribution
            self.direction_str = {'l2r', 'r2l'};
            self.next_ueId = 1; 
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

            max_mMTC_packet_size = Constants_c.max_mMTC_packet_size; 
            log2_max_mMTC_packet_size = log2(max_mMTC_packet_size);
            % Here is a simple solution, where it is to use the multinomial
            % distribution with equal probability, i.e. uniform sampling
            % between 1 and log2_max_mMTC_packet_size)
            self.pd_mMTC_log2_packet_size = makedist('Multinomial', 'Probabilities', ...
                1 / log2_max_mMTC_packet_size * ones(1, log2_max_mMTC_packet_size));  

            p_mMTC_tx = Constants_c.p_mMTC_tx;
            self.pd_mMTC_tx = makedist('Binomial', 'N', 1, 'p', p_mMTC_tx);

            self.winner2_channel_smp_time = 1600; % recommeded value from mathworld's default configuration

            self.nw_characteristics_v2x = struct('SINRdB', zeros(self.sim_nr_drops_per_round, self.M_veh), 'nr_packets', zeros(self.sim_nr_drops_per_round, self.M_veh));
            self.nw_characteristics_eMBB = struct('SINRdB', zeros(self.sim_nr_drops_per_round, self.M_eMBB), 'nr_session', zeros(self.sim_nr_drops_per_round, self.M_eMBB));
            self.nw_characteristics_mMTC = struct('SINRdB', zeros(self.sim_nr_drops_per_round, self.M_mMTC), 'packet_size', zeros(self.sim_nr_drops_per_round, self.M_mMTC), 'tx_flags', zeros(self.sim_nr_drops_per_round, self.M_mMTC));

            self.prb_utilization_ratio_v2x = zeros(self.sim_nr_drops_per_round, 1);
            self.prb_utilization_ratio_eMBB = zeros(self.sim_nr_drops_per_round, 1);
            self.prb_utilization_ratio_mMTC = zeros(self.sim_nr_drops_per_round, 1);
            self.prb_utilization_ratio_total = zeros(self.sim_nr_drops_per_round, 1);
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

        function [] = init_slicing_ratio(self, alpha_v2x_i, alpha_eMBB_i, alpha_mMTC_i)
            self.alpha_v2x = alpha_v2x_i;
            self.alpha_eMBB = alpha_eMBB_i;
            self.alpha_mMTC = alpha_mMTC_i;
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
                y = self.gNB_Bs.pos_y + r * sin(theta);
                self.eMBB_UEs{count} = eMBB_UE_c(self.next_ueId, x, y, ...
                    Constants_c.UE_height);
                self.next_ueId = self.next_ueId + 1;
            end
            % mMTC UEs are uniformly distribution around the base station
            % with a minimum and maximum distrance
            self.mMTC_UEs = cell(self.M_mMTC, 1);
            for count = 1 : self.M_mMTC
                r = self.pd_radius.random();
                theta = self.pd_theta.random();
                x = self.gNB_Bs.pos_x + r * cos(theta);
                y = self.gNB_Bs.pos_y + r * sin(theta);
                log2_packet_size = self.pd_mMTC_log2_packet_size.random();
                packet_size = 2 ^ log2_packet_size;
                tx_flag = self.pd_mMTC_tx.random();
                self.mMTC_UEs{count} = mMTC_UE_c(self.next_ueId, x, y, ...
                    Constants_c.UE_height, packet_size, tx_flag);
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
            % rate
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

            % For each mMTC UE, it should update its wakeup-transmission
            % flag
            for count = 1 : self.M_mMTC
                self.mMTC_UEs{count}.tx_flag = self.pd_mMTC_tx.random();
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
                plot([left(1), right(1)], [left(2), right(2)], '-k');
            end
          
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
                plot(x, y, "--y");
                text(x(2),y(2), sprintf('sinr=%.1f dB', self.mMTC_UEs{count}.sinr_dB));
            end            
        end

        function [] = run_winner2_model_simulation_drop(self, t_i)
            % should be override by its subclass
        end

        function [] = collect_nw_characteristic(self, t_i)
            for count = 1 : self.M_veh
                self.nw_characteristics_v2x.SINRdB(t_i, count) = self.vehicular_UEs{count}.sinr_dB;
                self.nw_characteristics_v2x.nr_packets(t_i, count) = self.vehicular_UEs{count}.nr_packet;
            end
            for count = 1 : self.M_eMBB
                self.nw_characteristics_eMBB.SINRdB(t_i, count) = self.eMBB_UEs{count}.sinr_dB;
                self.nw_characteristics_eMBB.nr_session(t_i, count) = self.eMBB_UEs{count}.nr_session();
            end
            for count = 1 : self.M_mMTC
                self.nw_characteristics_mMTC.SINRdB(t_i, count) = self.mMTC_UEs{count}.sinr_dB;
                self.nw_characteristics_mMTC.packet_size(t_i, count) = self.mMTC_UEs{count}.packet_size;
                self.nw_characteristics_mMTC.tx_flags(t_i, count) = self.mMTC_UEs{count}.tx_flag;
            end
        end

        function [] = calculate_PRB_utilization(self, t_i, dl_or_ul_i)
            t = t_i;
            dl_or_ul = dl_or_ul_i;
            Fd = Constants_c.Fd; 
            Sm = Constants_c.Sm;    
            B = Constants_c.B;
            % calculate PRB utilization for v2x slice
            Gamma1 = 0;
            Gamma1_max = self.alpha_v2x * Constants_c.N_RB;
            for count = 1 : self.M_veh
                m = self.vehicular_UEs{count}.nr_packet;
                SINR_dB = self.vehicular_UEs{count}.sinr_dB;
                Speff = mapping_SINR2Speff_f(SINR_dB, dl_or_ul);
                 
                if (Speff > 0)
                    Gamma1 = Gamma1 + m * Sm / (Speff * B * Fd);
                end                
            end   
            nr_PRB_v2x = min(Gamma1, Gamma1_max);

            % calculate PRB utilization for eMBB slice
            Gamma2 = 0;
            Gamma2_max = self.alpha_eMBB * Constants_c.N_RB;
            for count = 1 : self.M_eMBB
                Rb = Constants_c.Rb_session * self.eMBB_UEs{count}.nr_session; 
                SINR_dB = self.eMBB_UEs{count}.sinr_dB;
                Speff = mapping_SINR2Speff_f(SINR_dB, dl_or_ul);
                if (Speff > 0)
                    rho = Rb / (Speff * B);
                    Gamma2 = Gamma2 + rho;
                end
            end
            nr_PRB_eMBB = min(Gamma2, Gamma2_max);

            % calculate PRB utilization for mMTC slice
            Gamma3 = 0;
            Gamma3_max = self.alpha_mMTC * Constants_c.N_RB;
            for count = 1 : self.M_mMTC
                tx = self.mMTC_UEs{count}.tx_flag;
                ps = self.mMTC_UEs{count}.packet_size;
                SINR_dB = self.mMTC_UEs{count}.sinr_dB;
                Speff = mapping_SINR2Speff_f(SINR_dB, dl_or_ul);
                if (Speff > 0)
                    Gamma3 = Gamma3 + tx * ps / (Speff * B * Fd);
                end
            end
            nr_PRB_mMTC = min(Gamma3, Gamma3_max);            

            self.prb_utilization_ratio_v2x(t) = nr_PRB_v2x / Gamma1_max;
            self.prb_utilization_ratio_eMBB(t) = nr_PRB_eMBB / Gamma2_max;
            self.prb_utilization_ratio_mMTC(t) = nr_PRB_mMTC / Gamma3_max;
            self.prb_utilization_ratio_total(t) = (nr_PRB_v2x + nr_PRB_eMBB + nr_PRB_mMTC) / (Gamma1_max + Gamma2_max + Gamma3_max);

            fprintf('PRB utilization at t=%d: %f/%f, %f/%f, %f/%f\n', t_i, nr_PRB_v2x, Gamma1_max, nr_PRB_eMBB, Gamma2_max, nr_PRB_mMTC, Gamma3_max);
        end

        function [avg_ratio_o] = calculate_avg_PRB_utilization(self)
            avg_ratio_o = mean(self.prb_utilization_ratio_total);
        end

        function [] = update_slicing_ratio(self, t_i)
            % SHOULD BE OVERRIDE BY ITS SUBCLASS
        end

        function [] = start_simulation(self)
            % SHOULD BE OVERRIDE BY ITS SUBCLASS
        end
    end

end