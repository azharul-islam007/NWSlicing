function [] = network_model_dl_demo_f()

Constants_c.print_configurations();
fprintf(1, 'start dl network model...\n');
M_eMBB = Constants_c.M_eMBB;
M_veh = Constants_c.M_veh;
M_mMTC = Constants_c.M_mMTC;

sim_nr_drops_per_round = Constants_c.sim_nr_drops_per_round;
sim_nr_round = Constants_c.sim_nr_round;

lambda_e_lst = Constants_c.lambda_e_lst;
avg_prb_utilization = zeros(length(lambda_e_lst), 1);
for count = 1 : length(lambda_e_lst)
    lambda_e = lambda_e_lst(count);
    % create the network model for simulation
    nw = Network_Model_dl_c(M_eMBB, M_veh, M_mMTC, lambda_e, sim_nr_drops_per_round, sim_nr_round);
    % start to run the network model
    avg_prb_utilization(count) = nw.start_simulation();
end

for count = 1 : length(lambda_e_lst)
    fprintf(1, '%f -> %f\n', lambda_e_lst(count), avg_prb_utilization(count));
end
end

