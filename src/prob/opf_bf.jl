""
function solve_opf_bf(file, model_type::Type{T}, optimizer; kwargs...) where T <: AbstractBFModel
    return solve_model(file, model_type, optimizer, build_opf_bf; kwargs...)
end

"""
Solve the OPF problem with proposed dynamic relaxation methods
Recommended arguments:
- `K_init` - initial value of level of R&f
- `K_max` - maximum value of level of R&f 
- `constraints_flag` - type of constraints
    - 0: line approximation
    - 1: region relaxation (basic)
    - 2: region relaxation + trivial cut (RRC)
    - 3: pyramidal relaxation (PR)
    - 4: dynamic pyramidal relaxation (DPR)
    - 5: MISOCP + DPR
    - 6: MISOCP + PR
- `outer_flag` - type of outer approximation (effective only for constraints_flag = 4,5)
    - 1: global outer approximation
    - 2: delayed outer approximation
    - 3: local outer approximation based on original variables
    - 4: local outer approximation based on initial folded variables
    - 5: local outer approximation based on latest folded variables
- `randseed` - random seed for the random number generator
- `PwdRatio` - ratio of wind power to the total power
"""
function solve_opf_bf_dr(file, model_type::Type{T}, optimizer; kwargs...) where T <: AbstractBFModel
    pm, result = get_and_solve_model(file, model_type, optimizer, build_opf_bf_dr; kwargs...)
    return pm, result
end

""
function solve_mn_opf_bf_strg(file, model_type::Type{T}, optimizer; kwargs...) where T <: AbstractBFModel
    return solve_model(file, model_type, optimizer, build_mn_opf_bf_strg; multinetwork=true, kwargs...)
end

""
function build_opf_bf(pm::AbstractPowerModel)
    variable_bus_voltage(pm)
    variable_gen_power(pm)
    variable_branch_power(pm)
    variable_branch_current(pm)
    variable_dcline_power(pm)

    # objective_min_fuel_and_flow_cost(pm)
    objective_min_linear_cost(pm)

    constraint_model_current(pm)

    for i in ids(pm, :ref_buses)
        constraint_theta_ref(pm, i)
    end

    for i in ids(pm, :bus)
        constraint_power_balance(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_power_losses(pm, i)
        constraint_voltage_magnitude_difference(pm, i)

        constraint_voltage_angle_difference(pm, i)

        constraint_thermal_limit_from(pm, i)
        # constraint_thermal_limit_to(pm, i)
    end

    for i in ids(pm, :dcline)
        constraint_dcline_power_losses(pm, i)
    end
end

"dr methods of the bf model"
function build_opf_bf_dr(pm::AbstractPowerModel)
    param = ref(pm, :param)
    
    K_init = haskey(param, "K_init") ? param["K_init"] : 0
    K_max = haskey(param, "K_max") ? param["K_max"] : 4
    constraints_flag = haskey(param, "constraints_flag") ? param["constraints_flag"] : 4
    outer_flag = haskey(param, "outer_flag") ? param["outer_flag"] : 5
    randseed = haskey(param, "randseed") ? param["randseed"] : 0
    PwdRatio = haskey(param, "PwdRatio") ? param["PwdRatio"] : 0

    println("K_init: ", K_init, " K_max: ", K_max, " constraints_flag: ", constraints_flag,
     " outer_flag: ", outer_flag, " randseed: ", randseed, " PwdRatio: ", PwdRatio)

    # Set the model buspair_parameters
    model = pm.model
    # JuMP.set_attribute(model, "Threads", 8)
    JuMP.set_attribute(model, "Method", 2)
    JuMP.set_attribute(model, "MIPGap", 1e-3)
    # JuMP.set_attribute(model, "TimeLimit", 7200)
    JuMP.set_attribute(model, "TimeLimit", 3600)      #temporary 
    JuMP.set_attribute(model, "Seed",randseed)
    # JuMP.set_attribute(model, "SolutionLimit",1)

    variable_bus_voltage(pm)
    variable_gen_power(pm)
    variable_branch_power(pm)
    variable_branch_current(pm)
    variable_dcline_power(pm)
    variable_auxiliary(pm, K_max)   # prepare the auxiliary variables for dynamic relaxation methods

    # objective_min_fuel_and_flow_cost(pm)
    objective_min_linear_cost(pm)

    # Replace the conic constraints with piecewise linear constraints
    # constraint_model_current(pm)
    constraint_axis_symmetric(pm, K_init, K_max, constraints_flag)
    constraint_rotation_and_fold(pm, K_init, K_max, constraints_flag)


    for i in ids(pm, :ref_buses)
        constraint_theta_ref(pm, i)
    end

    for i in ids(pm, :bus)
        constraint_power_balance(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_power_losses(pm, i)
        constraint_voltage_magnitude_difference(pm, i)

        constraint_voltage_angle_difference(pm, i)

        # Replace the conic constraints with piecewise linear constraints
        # constraint_thermal_limit_from(pm, i)
        # constraint_thermal_limit_to(pm, i)
    end

    for i in ids(pm, :dcline)
        constraint_dcline_power_losses(pm, i)
    end

    # Add the dynamic relaxation methods
    if constraints_flag == 4 || constraints_flag == 5
        callback_setup(pm, K_init, K_max, constraints_flag, outer_flag)
    end
end

"Build multinetwork branch flow storage OPF"
function build_mn_opf_bf_strg(pm::AbstractPowerModel)
    for (n, network) in nws(pm)
        variable_bus_voltage(pm, nw=n)
        variable_gen_power(pm, nw=n)
        variable_storage_power_mi(pm, nw=n)
        variable_branch_power(pm, nw=n)
        variable_branch_current(pm, nw=n)
        variable_dcline_power(pm, nw=n)

        constraint_model_current(pm, nw=n)

        for i in ids(pm, :ref_buses, nw=n)
            constraint_theta_ref(pm, i, nw=n)
        end

        for i in ids(pm, :bus, nw=n)
            constraint_power_balance(pm, i, nw=n)
        end

        for i in ids(pm, :storage, nw=n)
            constraint_storage_complementarity_mi(pm, i, nw=n)
            constraint_storage_losses(pm, i, nw=n)
            constraint_storage_thermal_limit(pm, i, nw=n)
        end

        for i in ids(pm, :branch, nw=n)
            constraint_power_losses(pm, i, nw=n)
            constraint_voltage_magnitude_difference(pm, i, nw=n)

            constraint_voltage_angle_difference(pm, i, nw=n)

            constraint_thermal_limit_from(pm, i, nw=n)
            constraint_thermal_limit_to(pm, i, nw=n)
        end

        for i in ids(pm, :dcline, nw=n)
            constraint_dcline_power_losses(pm, i, nw=n)
        end
    end

    network_ids = sort(collect(nw_ids(pm)))

    n_1 = network_ids[1]
    for i in ids(pm, :storage, nw=n_1)
        constraint_storage_state(pm, i, nw=n_1)
    end

    for n_2 in network_ids[2:end]
        for i in ids(pm, :storage, nw=n_2)
            constraint_storage_state(pm, i, n_1, n_2)
        end
        n_1 = n_2
    end

    objective_min_fuel_and_flow_cost(pm)
end
