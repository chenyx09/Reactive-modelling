Code for the paper "Counter-example Guided Learning of Bounds on Environment Behavior"

Platform: Matlab
Prerequisite: Cplex, S-Taliro, Gurobi, qpOASES

To reproduce the result for the two robot problem, run 'robot_falsification'
To reproduce the result for the lane change problem, run 'LC_falsification'

Index of important functions:

params.m:                                     Declaring parameters
SVM_reactive_bound.m:                         Main SVM function that generates the classifier
extract_reactive_bound_2region_soft_bdry.m:   Generate reactive bound function from the classifier
BlackBoxMPCLaneChange.m:                      Lane change simulator
two_robot_sim.m:                              Simulation of the two robot example to generate positive data
two_robot_staliro_shortest_time_reactive.m:   S-Taliro file to generate negative data for the two robot example
lane_change_staliro.m:                        S-Taliro file to generate negative data for the lane change example