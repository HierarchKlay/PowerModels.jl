# this file contains dynamic relaxation formulation of branch flow model
using JuMP
"""
Defines relationship between branch (series) power flow, branch (series) current and node voltage magnitude
"""
function constraint_model_current(pm::AbstractSOCDRBFModel, n::Int)
    _check_missing_keys(var(pm, n), [:p,:q,:w,:ccm], typeof(pm))

    p  = var(pm, n, :p)
    q  = var(pm, n, :q)
    w  = var(pm, n, :w)
    ccm = var(pm, n, :ccm)

    for (i,branch) in ref(pm, n, :branch)
        f_bus = branch["f_bus"]
        t_bus = branch["t_bus"]
        f_idx = (i, f_bus, t_bus)
        tm = branch["tap"]

        # JuMP.@constraint(pm.model, p[f_idx]^2 + q[f_idx]^2 <= (w[f_bus]/tm^2)*ccm[i])

    end
end

function variable_auxiliary(pm::AbstractSOCDRBFModel, K, nw::Int = nw_id_default, report::Bool = true)
    # apparent power flow of the branch
    s = var(pm, nw)[:s] = JuMP.@variable(pm.model,
        [(l,i,j) in ref(pm, nw, :arcs)], base_name="$(nw)_s", lower_bound=0.0
    )

    for i in ids(pm, :branch)
        branch = ref(pm, nw, :branch, i)
        f_bus = branch["f_bus"]
        t_bus = branch["t_bus"]
        f_idx = (i, f_bus, t_bus)
        t_idx = (i, t_bus, f_bus)

        if haskey(branch, "rate_a")
            # build upper bound for variable s 
            JuMP.set_upper_bound(s[f_idx], branch["rate_a"])
            JuMP.set_upper_bound(s[t_idx], branch["rate_a"])
        end
    end

    # rotated second order cone variables
    vpi = var(pm, nw)[:vpi] = JuMP.@variable(pm.model,
        [(l,i,j) in ref(pm, nw, :arcs)], base_name="$(nw)_vpi", lower_bound=0.0
    )
    vmi = var(pm, nw)[:vmi] = JuMP.@variable(pm.model,
        [(l,i,j) in ref(pm, nw, :arcs)], base_name="$(nw)_vmi"
    )

    # auxiliary variables for axis symmetry
    pp = var(pm, nw)[:pp] = JuMP.@variable(pm.model,
        [(l,i,j) in ref(pm, nw, :arcs)], base_name="$(nw)_pp", lower_bound=0.0
    ) 
    pn = var(pm, nw)[:pn] = JuMP.@variable(pm.model,
        [(l,i,j) in ref(pm, nw, :arcs)], base_name="$(nw)_pn", lower_bound=0.0
    ) 
    flag_p = var(pm, nw)[:flag_p] = JuMP.@variable(pm.model,
        [(l,i,j) in ref(pm, nw, :arcs)], base_name="$(nw)_flag_p", Bin
    )
    qp = var(pm, nw)[:qp] = JuMP.@variable(pm.model,
        [(l,i,j) in ref(pm, nw, :arcs)], base_name="$(nw)_qp", lower_bound=0.0
    )
    qn = var(pm, nw)[:qn] = JuMP.@variable(pm.model,
        [(l,i,j) in ref(pm, nw, :arcs)], base_name="$(nw)_qn", lower_bound=0.0
    )
    flag_q = var(pm, nw)[:flag_q] = JuMP.@variable(pm.model,
        [(l,i,j) in ref(pm, nw, :arcs)], base_name="$(nw)_flag_q", Bin
    )
    vmip = var(pm, nw)[:vmip] = JuMP.@variable(pm.model,
        [(l,i,j) in ref(pm, nw, :arcs)], base_name="$(nw)_vmip", lower_bound=0.0
    )
    vmin = var(pm, nw)[:vmin] = JuMP.@variable(pm.model,
        [(l,i,j) in ref(pm, nw, :arcs)], base_name="$(nw)_vmin", lower_bound=0.0
    )
    flag_vmi = var(pm, nw)[:flag_vmi] = JuMP.@variable(pm.model,
        [(l,i,j) in ref(pm, nw, :arcs)], base_name="$(nw)_flag_vmi", Bin
    )

    # auxiliary variables for R&F: P^2 +Q^2 <= S^2
    gs = []; hs = []; hs_temp = []; hs_temp_p = []; hs_temp_n = []; flag_hs_temp = []; 
    for k in 1:K+1
        gsk = JuMP.@variable(pm.model,
            [(l,i,j) in ref(pm, nw, :arcs)], base_name="$(nw)_gs_$(k)"
        )
        push!(gs, gsk)
        hsk = JuMP.@variable(pm.model,
            [(l,i,j) in ref(pm, nw, :arcs)], base_name="$(nw)_hs_$(k)"
        )
        push!(hs, hsk)
        hsk_temp = JuMP.@variable(pm.model,
            [(l,i,j) in ref(pm, nw, :arcs)], base_name="$(nw)_hs_temp_$(k)"
        )
        push!(hs_temp, hsk_temp)
        hsk_temp_p = JuMP.@variable(pm.model,
            [(l,i,j) in ref(pm, nw, :arcs)], base_name="$(nw)_hs_temp_p_$(k)", lower_bound=0.0
        )
        push!(hs_temp_p, hsk_temp_p)
        hsk_temp_n = JuMP.@variable(pm.model,
            [(l,i,j) in ref(pm, nw, :arcs)], base_name="$(nw)_hs_temp_n_$(k)", lower_bound=0.0
        )
        push!(hs_temp_n, hsk_temp_n)
        flag_hsk_temp = JuMP.@variable(pm.model,
            [(l,i,j) in ref(pm, nw, :arcs)], base_name="$(nw)_flag_hs_temp_$(k)", Bin
        )
        push!(flag_hs_temp, flag_hsk_temp)
    end

    # auxiliary variables for R&F: vmi^2 + S^2 <= vpi^2
    gi = []; hi = []; hi_temp = []; hi_temp_p = []; hi_temp_n = []; flag_hi_temp = [];
    for k in 1:K+1
        gik = JuMP.@variable(pm.model,
            [(l,i,j) in ref(pm, nw, :arcs)], base_name="$(nw)_gi_$(k)"
        )
        push!(gi, gik)
        hik = JuMP.@variable(pm.model,
            [(l,i,j) in ref(pm, nw, :arcs)], base_name="$(nw)_hi_$(k)"
        )
        push!(hi, hik)
        hik_temp = JuMP.@variable(pm.model,
            [(l,i,j) in ref(pm, nw, :arcs)], base_name="$(nw)_hi_temp_$(k)"
        )
        push!(hi_temp, hik_temp)
        hik_temp_p = JuMP.@variable(pm.model,
            [(l,i,j) in ref(pm, nw, :arcs)], base_name="$(nw)_hi_temp_p_$(k)", lower_bound=0.0
        )
        push!(hi_temp_p, hik_temp_p)
        hik_temp_n = JuMP.@variable(pm.model,
            [(l,i,j) in ref(pm, nw, :arcs)], base_name="$(nw)_hi_temp_n_$(k)", lower_bound=0.0
        )
        push!(hi_temp_n, hik_temp_n)
        flag_hik_temp = JuMP.@variable(pm.model,
            [(l,i,j) in ref(pm, nw, :arcs)], base_name="$(nw)_flag_hi_temp_$(k)", Bin
        )
        push!(flag_hi_temp, flag_hik_temp)
    end
    
    report && sol_component_value_edge(pm, nw, :branch, :sr, :st, ref(pm, nw, :arcs_from), ref(pm, nw, :arcs_to), s)
    pm.model[:vpi] = vpi; pm.model[:vmi] = vmi; 
    pm.model[:pp] = pp; pm.model[:pn] = pn; pm.model[:flag_p] = flag_p;
    pm.model[:qp] = qp; pm.model[:qn] = qn; pm.model[:flag_q] = flag_q;
    pm.model[:vmip] = vmip; pm.model[:vmin] = vmin; pm.model[:flag_vmi] = flag_vmi;
    pm.model[:gs] = gs; pm.model[:hs] = hs; pm.model[:hs_temp] = hs_temp;
    pm.model[:hs_temp_p] = hs_temp_p; pm.model[:hs_temp_n] = hs_temp_n; pm.model[:flag_hs_temp] = flag_hs_temp;
    pm.model[:gi] = gi; pm.model[:hi] = hi; pm.model[:hi_temp] = hi_temp;
    pm.model[:hi_temp_p] = hi_temp_p; pm.model[:hi_temp_n] = hi_temp_n; pm.model[:flag_hi_temp] = flag_hi_temp;
    
end


"define axis symmetric constraints for R&F"
function constraint_axis_symmetric(pm::AbstractSOCDRBFModel, K_init, K_max, constraints_flag, n::Int = nw_id_default, report::Bool = true)
    # check missing keys
    # _check_missing_keys(var(pm, nw), [:p,:q,:w,:ccm], typeof(pm))
    p = var(pm, n, :p); q = var(pm, n, :q); w = var(pm, n, :w); ccm = var(pm, n, :ccm); s = var(pm, n, :s); vpi = var(pm, n, :vpi); vmi = var(pm, n, :vmi)

    pp = var(pm, n, :pp); pn = var(pm, n, :pn); flag_p = var(pm, n, :flag_p)
    qp = var(pm, n, :qp); qn = var(pm, n, :qn); flag_q = var(pm, n, :flag_q)
    vmip = var(pm, n, :vmip); vmin = var(pm, n, :vmin); flag_vmi = var(pm, n, :flag_vmi)

    gs = pm.model[:gs]; hs = pm.model[:hs]
    hs_temp = pm.model[:hs_temp]; hs_temp_p = pm.model[:hs_temp_p]; hs_temp_n = pm.model[:hs_temp_n]; flag_hs_temp = pm.model[:flag_hs_temp]

    gi = pm.model[:gi]; hi = pm.model[:hi]
    hi_temp = pm.model[:hi_temp]; hi_temp_p = pm.model[:hi_temp_p]; hi_temp_n = pm.model[:hi_temp_n]; flag_hi_temp = pm.model[:flag_hi_temp]


    for (i,branch) in ref(pm, n, :branch)
        f_bus = branch["f_bus"]
        t_bus = branch["t_bus"]
        f_idx = (i, f_bus, t_bus)
        tm = branch["tap"]

        # axis symmetric constraints for R&F: P^2 + Q^2 <= S^2
        JuMP.@constraint(pm.model, p[f_idx] == pp[f_idx] - pn[f_idx], base_name = "p_decomp_$(i)")
        JuMP.@constraint(pm.model, gs[1][f_idx] == pp[f_idx] + pn[f_idx], base_name = "p_abs_$(i)")
        JuMP.@constraint(pm.model, q[f_idx] == qp[f_idx] - qn[f_idx], base_name = "q_decomp_$(i)")
        JuMP.@constraint(pm.model, hs[1][f_idx] == qp[f_idx] + qn[f_idx], base_name = "q_abs_$(i)")
        JuMP.@constraint(pm.model, pp[f_idx] <= branch["rate_a"] * flag_p[f_idx], base_name = "pp_bin_$(i)")
        JuMP.@constraint(pm.model, pn[f_idx] <= branch["rate_a"] * (1 - flag_p[f_idx]), base_name = "pn_bin_$(i)")
        JuMP.@constraint(pm.model, qp[f_idx] <= branch["rate_a"] * flag_q[f_idx], base_name = "qp_bin_$(i)")
        JuMP.@constraint(pm.model, qn[f_idx] <= branch["rate_a"] * (1 - flag_q[f_idx]), base_name = "qn_bin_$(i)")

        # axis symmetric constraints for R&F: vmi^2 + S^2 <= vpi^2
        bus = ref(pm, n, :bus, f_bus)
        vmi_min = 1/2 * (bus["vmin"]^2 / tm^2 - (branch["rate_a"]^2 * tm^2)/ bus["vmin"]^2)
        vmi_max = bus["vmax"]^2 / (2 *tm^2)
        Mh0 = max(0, -vmi_min)
        hmax = max(vmi_max, -vmi_min)
        M2 = sqrt(hmax^2 + branch["rate_a"]^2)
        JuMP.@constraint(pm.model, vmi[f_idx] == (w[f_bus]/tm^2 - ccm[i])/2, base_name = "vmi_def_$(i)")
        JuMP.@constraint(pm.model, vpi[f_idx] == (w[f_bus]/tm^2 + ccm[i])/2, base_name = "vpi_def_$(i)")
        JuMP.@constraint(pm.model, gi[1][f_idx] ==  s[f_idx], base_name = "s_abs_$(i)")
        JuMP.@constraint(pm.model, vmi[f_idx] == vmip[f_idx] - vmin[f_idx], base_name = "vmi_decomp_$(i)")
        JuMP.@constraint(pm.model, hi[1][f_idx] == vmip[f_idx] + vmin[f_idx], base_name = "vmi_abs_$(i)")
        JuMP.@constraint(pm.model, vmip[f_idx] <= vmi_max * flag_vmi[f_idx], base_name = "vmip_bin_$(i)")
        JuMP.@constraint(pm.model, vmin[f_idx] <= Mh0 * (1 - flag_vmi[f_idx]), base_name = "vmin_bin_$(i)")

    end
end

"define (pre-added) rotation and fold constraints for R&F"
function constraint_rotation_and_fold(pm::AbstractSOCDRBFModel, K_init, K_max, constraints_flag, n::Int = nw_id_default, report::Bool = true)
    # check missing keys
    # _check_missing_keys(var(pm, nw), [:p,:q,:w,:ccm], typeof(pm))
    p = var(pm, n, :p); q = var(pm, n, :q); w = var(pm, n, :w); ccm = var(pm, n, :ccm); s = var(pm, n, :s); vpi = var(pm, n, :vpi); vmi = var(pm, n, :vmi)

    pp = var(pm, n, :pp); pn = var(pm, n, :pn); flag_p = var(pm, n, :flag_p)
    qp = var(pm, n, :qp); qn = var(pm, n, :qn); flag_q = var(pm, n, :flag_q)
    vmip = var(pm, n, :vmip); vmin = var(pm, n, :vmin); flag_vmi = var(pm, n, :flag_vmi)

    gs = pm.model[:gs]; hs = pm.model[:hs]
    hs_temp = pm.model[:hs_temp]; hs_temp_p = pm.model[:hs_temp_p]; hs_temp_n = pm.model[:hs_temp_n]; flag_hs_temp = pm.model[:flag_hs_temp]

    gi = pm.model[:gi]; hi = pm.model[:hi]
    hi_temp = pm.model[:hi_temp]; hi_temp_p = pm.model[:hi_temp_p]; hi_temp_n = pm.model[:hi_temp_n]; flag_hi_temp = pm.model[:flag_hi_temp]


    for (i,branch) in ref(pm, n, :branch)
        f_bus = branch["f_bus"]
        t_bus = branch["t_bus"]
        f_idx = (i, f_bus, t_bus)
        tm = branch["tap"]

        # rotation-and-fold constraints: P^2 + Q^2 <= S^2
        if constraints_flag == 4 || constraints_flag == 5 # R&F constraints setup for dynamic methods
            if K_init == -1 && K > 0 # although add no initial cuts, basic R&F constraints are needed due to K > 0
                # R&F step 1
                k = 1
                # rotation
                JuMP.@constraint(pm.model, gs[k+1][f_idx] == gs[k][f_idx] * cos(pi/2^(k+1)) + hs[k][f_idx] * sin(pi/2^(k+1)), base_name = "rot_gs_$(k)_$(i)")
                JuMP.@constraint(pm.model, hs_temp[k][f_idx] == -gs[k][f_idx] * sin(pi/2^(k+1)) + hs[k][f_idx] * cos(pi/2^(k+1)), base_name = "rot_hs_$(k)_$(i)")
                # fold: continuous part
                JuMP.@constraint(pm.model, hs_temp[k][f_idx] == hs_temp_p[k][f_idx] - hs_temp_n[k][f_idx], base_name = "hs_temp_decomp_$(k)_$(i)")
                JuMP.@constraint(pm.model, hs[k+1][f_idx] == hs_temp_p[k][f_idx] + hs_temp_n[k][f_idx], base_name = "hs_abs_$(k)_$(i)")
                # fold: binary part
                JuMP.@constraint(pm.model, hs_temp_p[k][f_idx] <= branch["rate_a"] * flag_hs_temp[k][f_idx], base_name = "hs_temp_p_bin_$(k)_$(i)")
                JuMP.@constraint(pm.model, hs_temp_n[k][f_idx] <= branch["rate_a"] * (1 - flag_hs_temp[k][f_idx]), base_name = "hs_temp_n_bin_$(k)_$(i)")
            else
                for k in 1:K_init
                    # rotation
                    JuMP.@constraint(pm.model, gs[k+1][f_idx] == gs[k][f_idx] * cos(pi/2^(k+1)) + hs[k][f_idx] * sin(pi/2^(k+1)), base_name = "rot_gs_$(k)_$(i)")
                    JuMP.@constraint(pm.model, hs_temp[k][f_idx] == -gs[k][f_idx] * sin(pi/2^(k+1)) + hs[k][f_idx] * cos(pi/2^(k+1)), base_name = "rot_hs_$(k)_$(i)")
                    # fold: continuous part
                    JuMP.@constraint(pm.model, hs_temp[k][f_idx] == hs_temp_p[k][f_idx] - hs_temp_n[k][f_idx], base_name = "hs_temp_decomp_$(k)_$(i)")
                    JuMP.@constraint(pm.model, hs[k+1][f_idx] == hs_temp_p[k][f_idx] + hs_temp_n[k][f_idx], base_name = "hs_abs_$(k)_$(i)")
                    # fold: binary part
                    JuMP.@constraint(pm.model, hs_temp_p[k][f_idx] <= branch["rate_a"] * flag_hs_temp[k][f_idx], base_name = "hs_temp_p_bin_$(k)_$(i)")
                    JuMP.@constraint(pm.model, hs_temp_n[k][f_idx] <= branch["rate_a"] * (1 - flag_hs_temp[k][f_idx]), base_name = "hs_temp_n_bin_$(k)_$(i)")
                end
            end
        else
            # all rotation-and-fold constraints are added for static methods
            for k in 1:K_max
                # rotation
                JuMP.@constraint(pm.model, gs[k+1][f_idx] == gs[k][f_idx] * cos(pi/2^(k+1)) + hs[k][f_idx] * sin(pi/2^(k+1)), base_name = "rot_gs_$(k)_$(i)")
                JuMP.@constraint(pm.model, hs_temp[k][f_idx] == -gs[k][f_idx] * sin(pi/2^(k+1)) + hs[k][f_idx] * cos(pi/2^(k+1)), base_name = "rot_hs_$(k)_$(i)")
                # fold: continuous part
                JuMP.@constraint(pm.model, hs_temp[k][f_idx] == hs_temp_p[k][f_idx] - hs_temp_n[k][f_idx], base_name = "hs_temp_decomp_$(k)_$(i)")
                JuMP.@constraint(pm.model, hs[k+1][f_idx] == hs_temp_p[k][f_idx] + hs_temp_n[k][f_idx], base_name = "hs_abs_$(k)_$(i)")
                # fold: binary part
                JuMP.@constraint(pm.model, hs_temp_p[k][f_idx] <= branch["rate_a"] * flag_hs_temp[k][f_idx], base_name = "hs_temp_p_bin_$(k)_$(i)")
                JuMP.@constraint(pm.model, hs_temp_n[k][f_idx] <= branch["rate_a"] * (1 - flag_hs_temp[k][f_idx]), base_name = "hs_temp_n_bin_$(k)_$(i)")
            end
        end
        # shape constraints for R&F: P^2 + Q^2 <= S^2
        if constraints_flag == 4 || constraints_flag == 5 # R&F constraints setup for dynamic methods
            if K_init == -1 # only add outer flag for K_init = -1
                JuMP.@constraint(pm.model, gs[1][f_idx] <= s[f_idx], base_name = "outer_gs_1_$(i)")
                JuMP.@constraint(pm.model, hs[1][f_idx] <= s[f_idx], base_name = "outer_hs_1_$(i)")
            else
                # add region cuts for K_init >= 0
                JuMP.@constraint(pm.model, s[f_idx] * cos(pi/2^(K_init+2)) <= gs[K_init+1][f_idx] * cos(pi/2^(K_init+2)) + hs[K_init+1][f_idx] * sin(pi/2^(K_init+2)), base_name = "inner_s_$(K_init+1)_$(i)")
                JuMP.@constraint(pm.model, gs[K_init+1][f_idx] <= s[f_idx], base_name = "outer_gs_$(K_init+1)_$(i)")
                JuMP.@constraint(pm.model, gs[K_init+1][f_idx] * cos(pi/2^(K_init+1)) + hs[K_init+1][f_idx] * sin(pi/2^(K_init+1)) <= s[f_idx], base_name = "outer_hs_$(K_init+1)_$(i)")
            end
        elseif constraints_flag == 3 || constraints_flag == 6
            # add region cuts for K_max
            JuMP.@constraint(pm.model, s[f_idx] * cos(pi/2^(K_max+2)) <= gs[K_max+1][f_idx] * cos(pi/2^(K_max+2)) + hs[K_max+1][f_idx] * sin(pi/2^(K_max+2)), base_name = "inner_$(K_max+1)_$(i)")
            JuMP.@constraint(pm.model, gs[K_max+1][f_idx] <= s[f_idx], base_name = "outer_gs_$(K_max+1)_$(i)")
            JuMP.@constraint(pm.model, gs[K_max+1][f_idx] * cos(pi/2^(K_max+1)) + hs[K_max+1][f_idx] * sin(pi/2^(K_max+1)) <= s[f_idx], base_name = "outer_hs_$(K_max+1)_$(i)")
        elseif constraints_flag == 2
            #TODO: add region cuts and trivial cuts for K_max
        elseif constraints_flag == 1
            #TODO: add region cuts for K_max
        elseif constraints_flag == 0
            JuMP.@constraint(pm.model, gs[K_max+1][f_idx] == s[f_idx] * cos(pi/2^(K_max+1)), base_name = "vertical_s_$(K_max+1)_$(i)")
            JuMP.@constraint(pm.model, hs[K_max+1][f_idx] >= 0, base_name = "positive_segment_s_$(K_max+1)_$(i)")
            JuMP.@constraint(pm.model, hs[K_max+1][f_idx] <= gs[K_max+1][f_idx] * tan(pi/2^(K_max+1)), base_name = "slope_s_$(K_max+1)_$(i)")
        end

        
        # axis symmetric constraints for R&F: vmi^2 + S^2 <= vpi^2
        bus = ref(pm, n, :bus, f_bus)
        vmi_min = 1/2 * (bus["vmin"]^2 / tm^2 - (branch["rate_a"]^2 * tm^2)/ bus["vmin"]^2)
        vmi_max = bus["vmax"]^2 / (2 * tm^2)
        Mh0 = max(0, -vmi_min)
        hmax = max(vmi_max, -vmi_min)
        M2 = sqrt(hmax^2 + branch["rate_a"]^2)
        # rotation-and-fold constraints: vmi^2 + S^2 <= vpi^2
        if constraints_flag == 4 || constraints_flag == 5 # R&F constraints setup for dynamic methods
            if K_init == -1 && K > 0 # although add no initial cuts, basic R&F constraints are needed due to K > 0
                # R&F step 1
                k = 1
                # rotation
                JuMP.@constraint(pm.model, gi[k+1][f_idx] == gi[k][f_idx] * cos(pi/2^(k+1)) + hi[k][f_idx] * sin(pi/2^(k+1)), base_name = "rot_gi_$(k)_$(i)")
                JuMP.@constraint(pm.model, hi_temp[k][f_idx] == -gi[k][f_idx] * sin(pi/2^(k+1)) + hi[k][f_idx] * cos(pi/2^(k+1)), base_name = "rot_hi_$(k)_$(i)")
                # fold: continuous part
                JuMP.@constraint(pm.model, hi_temp[k][f_idx] == hi_temp_p[k][f_idx] - hi_temp_n[k][f_idx], base_name = "hi_temp_decomp_$(k)_$(i)")
                JuMP.@constraint(pm.model, hi[k+1][f_idx] == hi_temp_p[k][f_idx] + hi_temp_n[k][f_idx], base_name = "hi_abs_$(k)_$(i)")
                # fold: binary part
                JuMP.@constraint(pm.model, hi_temp_p[k][f_idx] <= M2 * flag_hi_temp[k][f_idx], base_name = "hi_temp_p_bin_$(k)_$(i)")
                JuMP.@constraint(pm.model, hi_temp_n[k][f_idx] <= M2 * (1 - flag_hi_temp[k][f_idx]), base_name = "hi_temp_n_bin_$(k)_$(i)")
            else
                for k in 1:K_init
                    # rotation
                    JuMP.@constraint(pm.model, gi[k+1][f_idx] == gi[k][f_idx] * cos(pi/2^(k+1)) + hi[k][f_idx] * sin(pi/2^(k+1)), base_name = "rot_gi_$(k)_$(i)")
                    JuMP.@constraint(pm.model, hi_temp[k][f_idx] == -gi[k][f_idx] * sin(pi/2^(k+1)) + hi[k][f_idx] * cos(pi/2^(k+1)), base_name = "rot_hi_$(k)_$(i)")
                    # fold: continuous part
                    JuMP.@constraint(pm.model, hi_temp[k][f_idx] == hi_temp_p[k][f_idx] - hi_temp_n[k][f_idx], base_name = "hi_temp_decomp_$(k)_$(i)")
                    JuMP.@constraint(pm.model, hi[k+1][f_idx] == hi_temp_p[k][f_idx] + hi_temp_n[k][f_idx], base_name = "hi_abs_$(k)_$(i)")
                    # fold: binary part
                    JuMP.@constraint(pm.model, hi_temp_p[k][f_idx] <= M2 * flag_hi_temp[k][f_idx], base_name = "hi_temp_p_bin_$(k)_$(i)")
                    JuMP.@constraint(pm.model, hi_temp_n[k][f_idx] <= M2 * (1 - flag_hi_temp[k][f_idx]), base_name = "hi_temp_n_bin_$(k)_$(i)")
                end
            end
        else
            # all rotation-and-fold constraints are added for static methods
            for k in 1:K_max
                # rotation
                JuMP.@constraint(pm.model, gi[k+1][f_idx] == gi[k][f_idx] * cos(pi/2^(k+1)) + hi[k][f_idx] * sin(pi/2^(k+1)), base_name = "rot_gi_$(k)_$(i)")
                JuMP.@constraint(pm.model, hi_temp[k][f_idx] == -gi[k][f_idx] * sin(pi/2^(k+1)) + hi[k][f_idx] * cos(pi/2^(k+1)), base_name = "rot_hi_$(k)_$(i)")
                # fold: continuous part
                JuMP.@constraint(pm.model, hi_temp[k][f_idx] == hi_temp_p[k][f_idx] - hi_temp_n[k][f_idx], base_name = "hi_temp_decomp_$(k)_$(i)")
                JuMP.@constraint(pm.model, hi[k+1][f_idx] == hi_temp_p[k][f_idx] + hi_temp_n[k][f_idx], base_name = "hi_abs_$(k)_$(i)")
                # fold: binary part
                JuMP.@constraint(pm.model, hi_temp_p[k][f_idx] <= M2 * flag_hi_temp[k][f_idx], base_name = "hi_temp_p_bin_$(k)_$(i)")
                JuMP.@constraint(pm.model, hi_temp_n[k][f_idx] <= M2 * (1 - flag_hi_temp[k][f_idx]), base_name = "hi_temp_n_bin_$(k)_$(i)")
            end
        end
        # shape constraints for R&F: vmi^2 + S^2 <= vpi^2
        if constraints_flag == 4 || constraints_flag == 5 # R&F constraints setup for dynamic methods
            if K_init == -1 # only add outer flag for K_init = -1
                JuMP.@constraint(pm.model, gi[1][f_idx] <= vpi[f_idx], base_name = "outer_gi_1_$(i)")
                JuMP.@constraint(pm.model, hi[1][f_idx] <= vpi[f_idx], base_name = "outer_hi_1_$(i)")
            else
                # add region cuts for K_init >= 0
                JuMP.@constraint(pm.model, vpi[f_idx] * cos(pi/2^(K_init+2)) <= gi[K_init+1][f_idx] * cos(pi/2^(K_init+2)) + hi[K_init+1][f_idx] * sin(pi/2^(K_init+2)), base_name = "inner_i_$(K_init+1)_$(i)")
                JuMP.@constraint(pm.model, gi[K_init+1][f_idx] <= vpi[f_idx], base_name = "outer_gi_$(K_init+1)_$(i)")
                JuMP.@constraint(pm.model, gi[K_init+1][f_idx] * cos(pi/2^(K_init+1)) + hi[K_init+1][f_idx] * sin(pi/2^(K_init+1)) <= vpi[f_idx], base_name = "outer_hi_$(K_init+1)_$(i)")
            end
        elseif constraints_flag == 3 || constraints_flag == 6
            # add region cuts for K_max
            JuMP.@constraint(pm.model, vpi[f_idx] * cos(pi/2^(K_max+2)) <= gi[K_max+1][f_idx] * cos(pi/2^(K_max+2)) + hi[K_max+1][f_idx] * sin(pi/2^(K_max+2)), base_name = "inner_$(K_max+1)_$(i)")
            JuMP.@constraint(pm.model, gi[K_max+1][f_idx] <= vpi[f_idx], base_name = "outer_gi_$(K_max+1)_$(i)")
            JuMP.@constraint(pm.model, gi[K_max+1][f_idx] * cos(pi/2^(K_max+1)) + hi[K_max+1][f_idx] * sin(pi/2^(K_max+1)) <= vpi[f_idx], base_name = "outer_hi_$(K_max+1)_$(i)")
        elseif constraints_flag == 2
            #TODO: add region cuts and trivial cuts for K_max
        elseif constraints_flag == 1
            #TODO: add region cuts for K_max
        elseif constraints_flag == 0
            JuMP.@constraint(pm.model, gi[K_max+1][f_idx] == vpi[f_idx] * cos(pi/2^(K_max+1)), base_name = "vertical_i_$(K_max+1)_$(i)")
            JuMP.@constraint(pm.model, hi[K_max+1][f_idx] >= 0, base_name = "positive_segment_i_$(K_max+1)_$(i)")
            JuMP.@constraint(pm.model, hi[K_max+1][f_idx] <= gi[K_max+1][f_idx] * tan(pi/2^(K_max+1)), base_name = "slope_i_$(K_max+1)_$(i)")
        end

        if constraints_flag == 5 || constraints_flag == 6
            JuMP.@constraint(pm.model, p[f_idx]^2 + q[f_idx]^2 <= s[f_idx]^2, base_name = "socp_1_$(i)")
            JuMP.@constraint(pm.model, s[f_idx]^2 + vmi[f_idx]^2 <= vpi[f_idx]^2, base_name = "socp_2_$(i)")
        end
    end
end

"setup callback function for R&F"
function callback_setup(pm::AbstractSOCDRBFModel, K_init, K, constraints_flag, outer_flag, n::Int = nw_id_default, report::Bool = true)
    nCuts = pm.model[:nCuts] = [0]
    nOAIters = pm.model[:nOAIters] = [0]
    num_bra = length(ref(pm, n, :branch))
    println("num_bra = ", num_bra)
    if K_init == 0-1
        count_rf = pm.model[:count_rf] = fill(K_init+1, num_bra, 2)     # max num of r&f for each conic surface constraints till now
    else
        count_rf = pm.model[:count_rf] = fill(K_init, num_bra, 2)         # max num of r&f for each conic surface constraints till now
    end
    inrecord_rf = pm.model[:inrecord_rf] = Vector{Vector{Int}}()         # record history of added inner cut for each conic surface constraints
    for i in 1:2*num_bra
        if K_init == 0-1
            push!(inrecord_rf, [])
        else
            push!(inrecord_rf, [K_init])
        end
    end                                         # 1:num_bra vector:Sr Cons;num_bra+1:2num_bra:VrI_add Cons
    outrecord_rf = pm.model[:outrecord_rf] = Vector{Vector{Int}}()        # record history of added outer cut for each conic surface constraints
    for i in 1:2*num_bra
        push!(outrecord_rf, [K_init])
    end    
    outrecord_pq = pm.model[:outrecord_pq] = ones(2 * num_bra, 2^(K+2))      # out record for method 3
    for i in [1, 2^K+1, 2*2^K+1, 3*2^K+1]       # K_init=1-1 provide part of outer cuts
        outrecord_pq[:, i] .= 0
    end
    outrecord = pm.model[:outrecord] = ones(2 * num_bra, 2^K+1)      # out record for method 4
    for i in [1, 2^K+1]       # K_init=1-1 provide part of outer cuts
        outrecord[:, i] .= 0
    end
    outtensor = pm.model[:outtensor] = [[ones(2^(K-i)+1) for i in 0:K] for _ in 1:2*num_bra]   # out record for method 4 (outer_flag = 5)
    for t in 1:2*num_bra
        for i in [1, 2^K+1]       # K_init=1-1 provide part of outer cuts
            outtensor[t][1][i] = 0
        end
    end
    updateouttensor(outtensor,num_bra,K)
    function lazyCons(cb_data)
        isLazy = true
        callback_function(cb_data, pm, isLazy, K_init, K, constraints_flag, outer_flag, n)
    end
    MOI.set(
        pm.model,
        MOI.LazyConstraintCallback(),
        lazyCons,
    )
end

# forward propogation update for outtensor
function updateouttensor(outtensor,num_bra,K)
    for i in 1:2 * num_bra
        for t in 1:K
            for j in 1:2^(K-t)+1
                if outtensor[i][t][j] == 0 && outtensor[i][t][2^(K-t+1)+2-j] == 0
                    outtensor[i][t+1][2^(K-t)+2-j] = 0    #if corresponding 2 outer cuts in the last r&f already added, no need to add the current cut
                end
            end
        end
    end
end

"callback function for R&F"
function callback_function(cb_data, pm, isLazy, K_init, K, constraints_flag, outer_flag, n)
    eps = 1e-6
    Cons = isLazy ? MOI.LazyConstraint : MOI.UserCut

    model = pm.model
    num_bra = length(ref(pm, n, :branch))
    nCuts = model[:nCuts]
    nOAIters = model[:nOAIters]
    count_rf = model[:count_rf]
    inrecord_rf = model[:inrecord_rf]
    outrecord_rf = model[:outrecord_rf]
    outrecord_pq = model[:outrecord_pq]
    outrecord = model[:outrecord]
    outtensor = model[:outtensor]

    # get the status of the current solution
    status = callback_node_status(cb_data, model)
    if status == MOI.CALLBACK_NODE_STATUS_INTEGER
        println("-----------callback----------------")
        # get varibales reference
        p = var(pm, n, :p); q = var(pm, n, :q); w = var(pm, n, :w); ccm = var(pm, n, :ccm); s = var(pm, n, :s); vpi = var(pm, n, :vpi); vmi = var(pm, n, :vmi)

        pp = var(pm, n, :pp); pn = var(pm, n, :pn); flag_p = var(pm, n, :flag_p)
        qp = var(pm, n, :qp); qn = var(pm, n, :qn); flag_q = var(pm, n, :flag_q)
        vmip = var(pm, n, :vmip); vmin = var(pm, n, :vmin); flag_vmi = var(pm, n, :flag_vmi)
    
        gs = pm.model[:gs]; hs = pm.model[:hs]
        hs_temp = pm.model[:hs_temp]; hs_temp_p = pm.model[:hs_temp_p]; hs_temp_n = pm.model[:hs_temp_n]; flag_hs_temp = pm.model[:flag_hs_temp]
    
        gi = pm.model[:gi]; hi = pm.model[:hi]
        hi_temp = pm.model[:hi_temp]; hi_temp_p = pm.model[:hi_temp_p]; hi_temp_n = pm.model[:hi_temp_n]; flag_hi_temp = pm.model[:flag_hi_temp]
        
        # get the values of the variables
        vals_gs = []; vals_hs = []; vals_gi = []; vals_hi = []
        for k in 1:K+1
            push!(vals_gs, callback_value.(cb_data, gs[k]))
            push!(vals_hs, callback_value.(cb_data, hs[k]))
            push!(vals_gi, callback_value.(cb_data, gi[k]))
            push!(vals_hi, callback_value.(cb_data, hi[k]))
        end
        vals_s = callback_value.(cb_data, s)
        vals_vpi = callback_value.(cb_data, vpi)

        # copy the values of the variables for calculation in the loop
        temp_vals_gs = deepcopy(vals_gs)
        temp_vals_hs = deepcopy(vals_hs)
        temp_vals_gi = deepcopy(vals_gi)
        temp_vals_hi = deepcopy(vals_hi)

        ###--- cone1 ---###
        # inner cut1
        for (i,branch) in ref(pm, n, :branch)
            f_bus = branch["f_bus"]
            t_bus = branch["t_bus"]
            f_idx = (i, f_bus, t_bus)
            tm = branch["tap"]

            for k in K_init+1:K
                # calculate temp_vals[k] based on temp_vals[k-1]
                if k > count_rf[i,1]
                    temp_vals_gs[k+1][f_idx] = temp_vals_gs[k][f_idx] * cos(pi/2^(k+1)) + temp_vals_hs[k][f_idx] * sin(pi/2^(k+1))
                    temp_vals_hs[k+1][f_idx] = abs(-temp_vals_gs[k][f_idx] * sin(pi/2^(k+1)) + temp_vals_hs[k][f_idx] * cos(pi/2^(k+1)))
                end

                if temp_vals_gs[k+1][f_idx]^2 + temp_vals_hs[k+1][f_idx]^2 + eps < vals_s[f_idx]^2
                    if !isempty(inrecord_rf[i])
                        if inrecord_rf[i][end] < K && count_rf[i,1] == K && k == K_init
                            if vals_s[f_idx] * cos(pi/2^(k+2)) > eps + temp_vals_gs[k+1][f_idx] * cos(pi/2^(k+2)) + temp_vals_hs[k+1][f_idx] * sin(pi/2^(k+2))
                                inner_cut = @build_constraint(
                                    s[f_idx] * cos(pi/2^(k+2)) <= gs[k+1][f_idx] * cos(pi/2^(k+2)) + hs[k+1][f_idx] * sin(pi/2^(k+2))
                                )
                                MOI.submit(
                                    model,
                                    Cons(cb_data),
                                    inner_cut
                                )

                                append!(inrecord_rf[i], k)  # record the history of adding inner cut 
                                break   # stop the inner r&f for the current conic surface constraint
                            end
                        end

                        if count_rf[i,1] < k    # check if the k-r&f mapping constraints are added
                            for tmp_k in count_rf[i,1]+1:k
                                # add r&f mapping constraints
                                cons_g = @build_constraint(gs[tmp_k+1][f_idx] == gs[tmp_k][f_idx] * cos(pi/2^(tmp_k+1)) + hs[tmp_k][f_idx] * sin(pi/2^(tmp_k+1)))
                                cons_h1 = @build_constraint(hs_temp[tmp_k][f_idx] == -gs[tmp_k][f_idx] * sin(pi/2^(tmp_k+1)) + hs[tmp_k][f_idx] * cos(pi/2^(tmp_k+1)))
                                cons_h2 = @build_constraint(hs_temp[tmp_k][f_idx] == hs_temp_p[tmp_k][f_idx] - hs_temp_n[tmp_k][f_idx])
                                cons_h3 = @build_constraint(hs[tmp_k+1][f_idx] == hs_temp_p[tmp_k][f_idx] + hs_temp_n[tmp_k][f_idx])
                                cons_h4 = @build_constraint(hs_temp_p[tmp_k][f_idx] <= branch["rate_a"] * flag_hs_temp[tmp_k][f_idx])
                                cons_h5 = @build_constraint(hs_temp_n[tmp_k][f_idx] <= branch["rate_a"] * (1 - flag_hs_temp[tmp_k][f_idx]))
                                cons_bound = @build_constraint(gs[tmp_k+1][f_idx] <= s[f_idx])     # this cut is the equivalent cut of the final cuts

                                MOI.submit(model, Cons(cb_data), cons_g)
                                MOI.submit(model, Cons(cb_data), cons_h1)
                                MOI.submit(model, Cons(cb_data), cons_h2)
                                MOI.submit(model, Cons(cb_data), cons_h3)
                                MOI.submit(model, Cons(cb_data), cons_h4)
                                MOI.submit(model, Cons(cb_data), cons_h5)


                                if outer_flag == 5
                                    if outtensor[i][tmp_k+1][1] == 1
                                        # outer cut method 5 needs r&f add the bound outer cut in the same time
                                        MOI.submit(model, Cons(cb_data), cons_bound)    
                                        outtensor[i][tmp_k+1][1] = 0
                                        updateouttensor(outtensor,num_bra,K)
                                    end
                                end

                                count_rf[i,1] = tmp_k   # update the max num of r&f 
                            end

                            # check if temp_vals[k+1] violate the inner cut. if violated, add the cut and break
                            if inrecord_rf[i][end] < k
                                if vals_s[f_idx] * cos(pi/2^(k+2)) > eps + temp_vals_gs[k+1][f_idx] * cos(pi/2^(k+2)) + temp_vals_hs[k+1][f_idx] * sin(pi/2^(k+2))
                                    inner_cut = @build_constraint(
                                        s[f_idx] * cos(pi/2^(k+2)) <= gs[k+1][f_idx] * cos(pi/2^(k+2)) + hs[k+1][f_idx] * sin(pi/2^(k+2))
                                    )
                                    MOI.submit(
                                        model,
                                        Cons(cb_data),
                                        inner_cut
                                    )

                                    append!(inrecord_rf[i], k)  # record the history of adding inner cut 
                                    break   # stop the inner r&f for the current conic surface constraint
                                end
                            end
                        end
                    else
                        if k > 0
                            for tmp_k in count_rf[i,1]+1:k
                                # add r&f mapping constraints
                                cons_g = @build_constraint(gs[tmp_k+1][f_idx] == gs[tmp_k][f_idx] * cos(pi/2^(tmp_k+1)) + hs[tmp_k][f_idx] * sin(pi/2^(tmp_k+1)))
                                cons_h1 = @build_constraint(hs_temp[tmp_k][f_idx] == -gs[tmp_k][f_idx] * sin(pi/2^(tmp_k+1)) + hs[tmp_k][f_idx] * cos(pi/2^(tmp_k+1)))
                                cons_h2 = @build_constraint(hs_temp[tmp_k][f_idx] == hs_temp_p[tmp_k][f_idx] - hs_temp_n[tmp_k][f_idx])
                                cons_h3 = @build_constraint(hs[tmp_k+1][f_idx] == hs_temp_p[tmp_k][f_idx] + hs_temp_n[tmp_k][f_idx])
                                cons_h4 = @build_constraint(hs_temp_p[tmp_k][f_idx] <= branch["rate_a"] * flag_hs_temp[tmp_k][f_idx])
                                cons_h5 = @build_constraint(hs_temp_n[tmp_k][f_idx] <= branch["rate_a"] * (1 - flag_hs_temp[tmp_k][f_idx]))
                                cons_bound = @build_constraint(gs[tmp_k+1][f_idx] <= s[f_idx])     # this cut is the equivalent cut of the final cuts

                                MOI.submit(model, Cons(cb_data), cons_g)
                                MOI.submit(model, Cons(cb_data), cons_h1)
                                MOI.submit(model, Cons(cb_data), cons_h2)
                                MOI.submit(model, Cons(cb_data), cons_h3)
                                MOI.submit(model, Cons(cb_data), cons_h4)
                                MOI.submit(model, Cons(cb_data), cons_h5)


                                if outer_flag == 5
                                    if outtensor[i][tmp_k+1][1] == 1
                                        # outer cut method 5 needs r&f add the bound outer cut in the same time
                                        MOI.submit(model, Cons(cb_data), cons_bound)    
                                        outtensor[i][tmp_k+1][1] = 0
                                        updateouttensor(outtensor,num_bra,K)
                                    end
                                end

                                count_rf[i,1] = tmp_k   # update the max num of r&f 
                            end
                        end
                        # check if 0-inner cuts are violated? if so, add the inner cuts
                        if vals_s[f_idx] * cos(pi/2^(k+2)) > eps + temp_vals_gs[k+1][f_idx] * cos(pi/2^(k+2)) + temp_vals_hs[k+1][f_idx] * sin(pi/2^(k+2))
                            inner_cut = @build_constraint(
                                s[f_idx] * cos(pi/2^(k+2)) <= gs[k+1][f_idx] * cos(pi/2^(k+2)) + hs[k+1][f_idx] * sin(pi/2^(k+2))
                            )
                            MOI.submit(
                                model,
                                Cons(cb_data),
                                inner_cut
                            )

                            append!(inrecord_rf[i], k)  # record the history of adding inner cut 
                            break   # stop the inner r&f for the current conic surface constraint
                        end

                    end
                end



                             

            end



        end
           
        # outer cut1
        if constraints_flag == 4
            if outer_flag == 1
                # outer cut1 ::method 1: add outer cut in r&f way
                for (i,branch) in ref(pm, n, :branch)
                    f_bus = branch["f_bus"]
                    t_bus = branch["t_bus"]
                    f_idx = (i, f_bus, t_bus)
                    tm = branch["tap"]

                    for k in K_init:K
                        if k > count_rf[i,1]
                            # calculate temp_vals[k] based on temp_vals[k-1]
                            temp_vals_gs[k+1][f_idx] = temp_vals_gs[k][f_idx] * cos(pi/2^(k+1)) + temp_vals_hs[k][f_idx] * sin(pi/2^(k+1))
                            temp_vals_hs[k+1][f_idx] = abs(-temp_vals_gs[k][f_idx] * sin(pi/2^(k+1)) + temp_vals_hs[k][f_idx] * cos(pi/2^(k+1)))
                        end

                        if k == count_rf[i,1] # need not to add r&f mapping and may need add outer cut 
                            if k == K_init || k == outrecord_rf[i][end]
                                continue        # k-outer cut already added
                            elseif (temp_vals_gs[k+1][f_idx] > vals_s[f_idx] + eps) ||
                                # check if point violate outer cuts 
                                (temp_vals_gs[k+1][f_idx] * cos(pi/2^(k+1)) + temp_vals_hs[k+1][f_idx] * sin(pi/2^(k+1)) > vals_s[f_idx] + eps)
                                
                                # k-outer cut is violated
                                outer_cut1 = @build_constraint(gs[k+1][f_idx] <= s[f_idx])
                                outer_cut2 = @build_constraint(
                                    gs[k+1][f_idx] * cos(pi/2^(k+1)) + hs[k+1][f_idx] * sin(pi/2^(k+1)) <= s[f_idx]
                                )

                                MOI.submit(model, Cons(cb_data), outer_cut1)
                                MOI.submit(model, Cons(cb_data), outer_cut2)

                                append!(outrecord_rf[i], k)  # record the history of adding outer cut
                                break
                            else
                                continue # k-outer cut satisfied, thus need not to be added
                            end
                        end

                        if temp_vals_gs[k+1][f_idx]^2 + temp_vals_hs[k+1][f_idx]^2 > vals_s[f_idx]^2 +eps  # check if point is located outside the cone 
                            if outrecord_rf[i][end] < K && count_rf[i,1] == K && k == K     # all r&f mappings are added but the outer cut is not added
                                if (temp_vals_gs[k+1][f_idx] > vals_s[f_idx] + eps) ||
                                    # check if point violate outer cuts 
                                    (temp_vals_gs[k+1][f_idx] * cos(pi/2^(k+1)) + temp_vals_hs[k+1][f_idx] * sin(pi/2^(k+1)) > vals_s[f_idx] + eps)
                                    
                                    # k-outer cut is violated
                                    outer_cut1 = @build_constraint(gs[k+1][f_idx] <= s[f_idx])
                                    outer_cut2 = @build_constraint(
                                        gs[k+1][f_idx] * cos(pi/2^(k+1)) + hs[k+1][f_idx] * sin(pi/2^(k+1)) <= s[f_idx]
                                    )

                                    MOI.submit(model, Cons(cb_data), outer_cut1)
                                    MOI.submit(model, Cons(cb_data), outer_cut2)

                                    append!(outrecord_rf[i], k)  # record the history of adding outer cut
                                    break
                                end
                            end

                            if count_rf[i,1] < k    # check if the k-r&f mapping constraints are added
                                for tmp_k in count_rf[i,1]+1:k
                                    # add r&f mapping constraints
                                    cons_g = @build_constraint(gs[tmp_k+1][f_idx] == gs[tmp_k][f_idx] * cos(pi/2^(tmp_k+1)) + hs[tmp_k][f_idx] * sin(pi/2^(tmp_k+1)))
                                    cons_h1 = @build_constraint(hs_temp[tmp_k][f_idx] == -gs[tmp_k][f_idx] * sin(pi/2^(tmp_k+1)) + hs[tmp_k][f_idx] * cos(pi/2^(tmp_k+1)))
                                    cons_h2 = @build_constraint(hs_temp[tmp_k][f_idx] == hs_temp_p[tmp_k][f_idx] - hs_temp_n[tmp_k][f_idx])
                                    cons_h3 = @build_constraint(hs[tmp_k+1][f_idx] == hs_temp_p[tmp_k][f_idx] + hs_temp_n[tmp_k][f_idx])
                                    cons_h4 = @build_constraint(hs_temp_p[tmp_k][f_idx] <= branch["rate_a"] * flag_hs_temp[tmp_k][f_idx])
                                    cons_h5 = @build_constraint(hs_temp_n[tmp_k][f_idx] <= branch["rate_a"] * (1 - flag_hs_temp[tmp_k][f_idx]))

                                    MOI.submit(model, Cons(cb_data), cons_g)
                                    MOI.submit(model, Cons(cb_data), cons_h1)
                                    MOI.submit(model, Cons(cb_data), cons_h2)
                                    MOI.submit(model, Cons(cb_data), cons_h3)
                                    MOI.submit(model, Cons(cb_data), cons_h4)
                                    MOI.submit(model, Cons(cb_data), cons_h5)

                                    count_rf[i,1] = tmp_k   # update the max num of r&f
                                end

                                # check if temp_vals[k+1] violate the outer cut. if violated, add the cut and break
                                if (temp_vals_gs[k+1][f_idx] > vals_s[f_idx] + eps) ||
                                    # check if point violate outer cuts 
                                    (temp_vals_gs[k+1][f_idx] * cos(pi/2^(k+1)) + temp_vals_hs[k+1][f_idx] * sin(pi/2^(k+1)) > vals_s[f_idx] + eps)
                                    
                                    # k-outer cut is violated
                                    outer_cut1 = @build_constraint(gs[k+1][f_idx] <= s[f_idx])
                                    outer_cut2 = @build_constraint(
                                        gs[k+1][f_idx] * cos(pi/2^(k+1)) + hs[k+1][f_idx] * sin(pi/2^(k+1)) <= s[f_idx]
                                    )

                                    MOI.submit(model, Cons(cb_data), outer_cut1)
                                    MOI.submit(model, Cons(cb_data), outer_cut2)

                                    append!(outrecord_rf[i], k)  # record the history of adding outer cut
                                    break
                                end
                            end
                        end
                    end
                end
            elseif outer_flag == 2
                # outer cut1 ::method 2: add central cut when violation occurs
                # calculate delta
                ref_delta = zeros(K+3)
                for k in 1:K+3
                    ref_delta[k]=(1/(cos(pi/2^(k+2))).^2)-1  
                end
                for (i,branch) in ref(pm, n, :branch)
                    f_bus = branch["f_bus"]
                    t_bus = branch["t_bus"]
                    f_idx = (i, f_bus, t_bus)
                    tm = branch["tap"]

                    for k in count_rf[i,1]:K
                        if outrecord_rf[i][end] >= K    # number of outer cuts reached limit
                            break
                        end
                        # calculate temp_vals[k] based on temp_vals[k-1]
                        if k > count_rf[i,1]
                            temp_vals_gs[k+1][f_idx] = temp_vals_gs[k][f_idx] * cos(pi/2^(k+1)) + temp_vals_hs[k][f_idx] * sin(pi/2^(k+1))
                            temp_vals_hs[k+1][f_idx] = abs(-temp_vals_gs[k][f_idx] * sin(pi/2^(k+1)) + temp_vals_hs[k][f_idx] * cos(pi/2^(k+1)))
                        end

                        if temp_vals_gs[k+1][f_idx]^2 + temp_vals_hs[k+1][f_idx]^2 > vals_s[f_idx]^2 + eps  # check if point is located outside the cone 
                           if k > count_rf[i,1]
                                # add r&f mapping constraints
                                cons_g = @build_constraint(gs[k+1][f_idx] == gs[k][f_idx] * cos(pi/2^(k+1)) + hs[k][f_idx] * sin(pi/2^(k+1)))
                                cons_h1 = @build_constraint(hs_temp[k][f_idx] == -gs[k][f_idx] * sin(pi/2^(k+1)) + hs[k][f_idx] * cos(pi/2^(k+1)))
                                cons_h2 = @build_constraint(hs_temp[k][f_idx] == hs_temp_p[k][f_idx] - hs_temp_n[k][f_idx])
                                cons_h3 = @build_constraint(hs[k+1][f_idx] == hs_temp_p[k][f_idx] + hs_temp_n[k][f_idx])
                                cons_h4 = @build_constraint(hs_temp_p[k][f_idx] <= branch["rate_a"] * flag_hs_temp[k][f_idx])
                                cons_h5 = @build_constraint(hs_temp_n[k][f_idx] <= branch["rate_a"] * (1 - flag_hs_temp[k][f_idx]))
                                cons_bound = @build_constraint(gs[k+1][f_idx] <= s[f_idx])     # this cut is the equivalent cut of the final cuts

                                MOI.submit(model, Cons(cb_data), cons_g)
                                MOI.submit(model, Cons(cb_data), cons_h1)
                                MOI.submit(model, Cons(cb_data), cons_h2)
                                MOI.submit(model, Cons(cb_data), cons_h3)
                                MOI.submit(model, Cons(cb_data), cons_h4)
                                MOI.submit(model, Cons(cb_data), cons_h5)
                                MOI.submit(model, Cons(cb_data), cons_bound)

                                count_rf[i,1] = k   # update the max num of r&f
                                append!(outrecord_rf[i], k)  # record the history of adding outer cut
                           end

                           if k == count_rf[i,1]
                                tmp_delta = (temp_vals_gs[k+1][f_idx]^2 + temp_vals_hs[k+1][f_idx]^2 - vals_s[f_idx]^2) / vals_s[f_idx]^2 # error measure

                                # prepare outer cuts
                                cons_tan1 = @build_constraint(gs[k+1][f_idx] * cos(pi/2^(k+2)) + hs[k+1][f_idx] * sin(pi/2^(k+2)) <= s[f_idx])
                                cons_tan2 = @build_constraint(gs[k+1][f_idx] * cos(pi/2^(k+3)) + hs[k+1][f_idx] * sin(pi/2^(k+3)) <= s[f_idx])
                                cons_tan3 = @build_constraint(gs[k+1][f_idx] * cos(3*pi/2^(k+3)) + hs[k+1][f_idx] * sin(3*pi/2^(k+3)) <= s[f_idx])
                                cons_tan4 = @build_constraint(gs[k+1][f_idx] * cos(pi/2^(k+4)) + hs[k+1][f_idx] * sin(pi/2^(k+4)) <= s[f_idx])
                                cons_tan5 = @build_constraint(gs[k+1][f_idx] * cos(3*pi/2^(k+4)) + hs[k+1][f_idx] * sin(3*pi/2^(k+4)) <= s[f_idx])
                                cons_tan6 = @build_constraint(gs[k+1][f_idx] * cos(5*pi/2^(k+4)) + hs[k+1][f_idx] * sin(5*pi/2^(k+4)) <= s[f_idx])
                                cons_tan7 = @build_constraint(gs[k+1][f_idx] * cos(7*pi/2^(k+4)) + hs[k+1][f_idx] * sin(7*pi/2^(k+4)) <= s[f_idx])

                                if tmp_delta > ref_delta[k+1]+1e-10 # central cut should be add
                                    if k+1 > K
                                        break
                                    else
                                        MOI.submit(model, Cons(cb_data), cons_tan1)
                                        append!(outrecord_rf[i], k+1)  # record the history of adding outer cut
                                        break
                                    end
                                elseif tmp_delta > ref_delta[k+2]+1e-10 # 2 subcentral cuts should be add

                                    if k+2 > K
                                        break
                                    else
                                        MOI.submit(model, Cons(cb_data), cons_tan1)
                                        MOI.submit(model, Cons(cb_data), cons_tan2)
                                        MOI.submit(model, Cons(cb_data), cons_tan3)
                                        append!(outrecord_rf[i], k+2)  # record the history of adding outer cut
                                        break
                                    end
                                elseif tmp_delta > ref_delta[k+3]+1e-10 # 4 subcentral cuts should be add
                                    if k+3 > K
                                        break
                                    else
                                        MOI.submit(model, Cons(cb_data), cons_tan1)
                                        MOI.submit(model, Cons(cb_data), cons_tan2)
                                        MOI.submit(model, Cons(cb_data), cons_tan3)
                                        MOI.submit(model, Cons(cb_data), cons_tan4)
                                        MOI.submit(model, Cons(cb_data), cons_tan5)
                                        MOI.submit(model, Cons(cb_data), cons_tan6)
                                        MOI.submit(model, Cons(cb_data), cons_tan7)
                                        append!(outrecord_rf[i], k+3)  # record the history of adding outer cut
                                        break
                                    end
                                else
                                    if k >= K-3
                                        break
                                    else
                                        continue
                                    end
                                end
                            end
                        end
                    end
                end
            elseif outer_flag == 3
                ## outer cut1 ::method 3: add cuts based on arguments of (P,Q)
                # calculate delta
                ref_delta = zeros(Float64,K+3)
                for k in 1:K+3
                    ref_delta[k]=(1/(cos(pi/2^(k+2))).^2)-1  
                end

                vals_p = callback_value.(cb_data, p)
                vals_q = callback_value.(cb_data, q)

                for (i,branch) in ref(pm, n, :branch)
                    f_bus = branch["f_bus"]
                    t_bus = branch["t_bus"]
                    f_idx = (i, f_bus, t_bus)
                    tm = branch["tap"]

                    if vals_p[f_idx]^2 + vals_q[f_idx]^2 > eps + vals_s[f_idx]^2
                        delta = (vals_p[f_idx]^2 + vals_q[f_idx]^2 - vals_s[f_idx]^2) / vals_s[f_idx]^2 # error measure

                        # calculate the argument of (P,Q)
                        if vals_p[f_idx] == 0   #argument is pi/2 or 3pi/2
                            if vals_q[f_idx] > 0
                                arg = pi/2
                            elseif vals_q[f_idx] < 0
                                arg = 3pi/2
                            else 
                                @assert vals_s[f_idx] != 0 "P,Q = 0, S not equal to 0"
                            end
                        elseif vals_q[f_idx] == 0   #argument is pi or 0
                            if vals_p[f_idx] > 0
                                arg = 0
                            elseif vals_p[f_idx] < 0
                                arg = pi
                            else 
                                @assert vals_s[f_idx] != 0 "P,Q = 0, S not equal to 0"
                            end
                        else
                            if vals_p[f_idx] > 0 && vals_q[f_idx] > 0
                                arg = atan(vals_q[f_idx]/vals_p[f_idx])
                            elseif vals_p[f_idx] < 0 && vals_q[f_idx] > 0
                                arg = pi + atan(vals_q[f_idx]/vals_p[f_idx])
                            elseif vals_p[f_idx] < 0 && vals_q[f_idx] < 0
                                arg = pi + atan(vals_q[f_idx]/vals_p[f_idx])
                            elseif vals_p[f_idx] > 0 && vals_q[f_idx] < 0
                                arg = 2pi + atan(vals_q[f_idx]/vals_p[f_idx])
                            end
                        end

                        for tmp_k in 1:K
                            if delta > ref_delta[tmp_k] + eps   #solution violates the neighbored tmp_k-cuts
                                arg_coef = arg/(pi/2^(tmp_k+1))+1
                                arg_l = Int(round(floor(arg_coef)))
                                arg_u = Int(round(ceil(arg_coef)))

                                if arg_u == 2^(tmp_k+2)+1
                                    if arg_l < arg_u
                                        arg_u = 1
                                    elseif arg_l == arg_u
                                        arg_u = 1
                                        arg_l = 1
                                    end
                                end

                                idx_l = Int(round((arg_l - 1) * 2^(K-tmp_k)+1))        #number of the floor cut in all cuts 
                                idx_u = Int(round((arg_u - 1) * 2^(K-tmp_k)+1))        #number of the ceil cut in all cuts 

                                if arg_l < arg_u || (arg_u == 1 && arg_l != 1)                 #solution is not located at the border
                                    # @assert outrecord_pq[i,idx_l] == 1 || outrecord_pq[i,idx_u] == 1 "violated with cut added"

                                    if outrecord_pq[i,idx_l] == 1   #the floor cut is not added yet
                                        # add the outer cut
                                        cut = @build_constraint(
                                            p[f_idx] * cos((idx_l-1) * pi/2^(K+1)) + q[f_idx] * sin((idx_l-1) * pi/2^(K+1)) <= s[f_idx]
                                        )
                                        MOI.submit(model, Cons(cb_data), cut)
                                        outrecord_pq[i,idx_l] = 0
                                    end
                                    if outrecord_pq[i,idx_u] == 1   #the ceil cut is not added yet
                                        # add the outer cut
                                        cut = @build_constraint(
                                            p[f_idx] * cos((idx_u-1) * pi/2^(K+1)) + q[f_idx] * sin((idx_u-1) * pi/2^(K+1)) <= s[f_idx]
                                        )
                                        MOI.submit(model, Cons(cb_data), cut)
                                        outrecord_pq[i,idx_u] = 0
                                    end
                                elseif arg_l == arg_u || (arg_u == 1 && arg_l == 1)              #solution is located at the border
                                    # @assert outrecord_pq[i,idx_l] == 1 || outrecord_pq[i,idx_u] == 1 "violated with cut added"
                                    if outrecord_pq[i,idx_l] == 1
                                        cut = @build_constraint(
                                            p[f_idx] * cos((idx_l-1) * pi/2^(K+1)) + q[f_idx] * sin((idx_l-1) * pi/2^(K+1)) <= s[f_idx]
                                        )
                                        MOI.submit(model, Cons(cb_data), cut)
                                        outrecord_pq[i,idx_l] = 0
                                    end
                                end

                                break
                            end
                        end
                    end
                end
            elseif outer_flag == 4
                # outer cut1 ::method 4: add cuts based on arguments of (gs0,hs0)
                # calculate delta
                ref_delta = zeros(Float64,K)
                for k in 1:K
                    ref_delta[k]=(1/(cos(pi/2^(k+2))).^2)-1
                end

                for (i,branch) in ref(pm, n, :branch)
                    f_bus = branch["f_bus"]
                    t_bus = branch["t_bus"]
                    f_idx = (i, f_bus, t_bus)
                    tm = branch["tap"]

                    if vals_gs[1][f_idx]^2 + vals_hs[1][f_idx]^2 > eps + vals_s[f_idx]^2
                        delta = (vals_gs[1][f_idx]^2 + vals_hs[1][f_idx]^2 - vals_s[f_idx]^2) / vals_s[f_idx]^2 # error measure

                        # calculate the argument of (gs0,hs0)
                        if vals_gs[1][f_idx] == 0   #argument is pi/2 or 3pi/2
                            arg = pi/2
                            if vals_hs[1][f_idx] == 0
                                @assert vals_s[f_idx] != 0 "gs0,hs0 = 0, S not equal to 0"
                            end
                        elseif vals_hs[1][f_idx] == 0   #argument is pi or 0
                            arg = 0
                            if vals_gs[1][f_idx] == 0
                                @assert vals_s[f_idx] != 0 "gs0,hs0 = 0, S not equal to 0"
                            end
                        else
                            arg = atan(vals_hs[1][f_idx]/vals_gs[1][f_idx])
                        end

                        for tmp_k in 1:K
                            if delta > ref_delta[tmp_k] + eps
                                arg_coef = arg/(pi/2^(tmp_k+1))+1
                                arg_l = Int(round(floor(arg_coef)))     #number of the floor cut in 2^(tmp_k+2) cuts
                                arg_u = Int(round(ceil(arg_coef)))      #number of the ceil cut in 2^(tmp_k+2) cuts

                                idx_l = Int(round((arg_l - 1) * 2^(K-tmp_k)+1))        #number of the floor cut in all cuts
                                idx_u = Int(round((arg_u - 1) * 2^(K-tmp_k)+1))        #number of the ceil cut in all cuts

                                if arg_l < arg_u
                                    # @assert outrecord[i,idx_l] == 1 || outrecord[i,idx_u] == 1 "violated with cut added"

                                    if outrecord[i,idx_l] == 1      #the floor cut is not added yet
                                        # add the outer cut
                                        cut = @build_constraint(
                                            gs[1][f_idx] * cos((idx_l-1) * pi/2^(K+1)) + hs[1][f_idx] * sin((idx_l-1) * pi/2^(K+1)) <= s[f_idx]
                                        )
                                        MOI.submit(model, Cons(cb_data), cut)
                                        outrecord[i,idx_l] = 0
                                    end

                                    if outrecord[i,idx_u] == 1      #the ceil cut is not added yet
                                        # add the outer cut
                                        cut = @build_constraint(
                                            gs[1][f_idx] * cos((idx_u-1) * pi/2^(K+1)) + hs[1][f_idx] * sin((idx_u-1) * pi/2^(K+1)) <= s[f_idx]
                                        )
                                        MOI.submit(model, Cons(cb_data), cut)
                                        outrecord[i,idx_u] = 0
                                    end
                                elseif arg_l == arg_u
                                    # @assert outrecord[i,idx_l] == 1 || outrecord[i,idx_u] == 1 "violated with cut added"
                                    cut = @build_constraint(
                                        gs[1][f_idx] * cos((idx_l-1) * pi/2^(K+1)) + hs[1][f_idx] * sin((idx_l-1) * pi/2^(K+1)) <= s[f_idx]
                                    )
                                    MOI.submit(model, Cons(cb_data), cut)
                                    outrecord[i,idx_l] = 0
                                end

                                break
                            end
                        end
                    end
                end
            elseif outer_flag == 5
                # outer cut1 ::method 5: add cuts based on arguments of (gsk,hsk) where k = count_rf
                # calculate delta
                ref_delta = zeros(Float64,K)
                for k in 1:K
                    ref_delta[k]=(1/(cos(pi/2^(k+2))).^2)-1
                end

                for (i,branch) in ref(pm, n, :branch)
                    f_bus = branch["f_bus"]
                    t_bus = branch["t_bus"]
                    f_idx = (i, f_bus, t_bus)
                    tm = branch["tap"]

                    if vals_gs[count_rf[i,1]+1][f_idx]^2 + vals_hs[count_rf[i,1]+1][f_idx]^2 > eps + vals_s[f_idx]^2
                        delta = (vals_gs[count_rf[i,1]+1][f_idx]^2 + vals_hs[count_rf[i,1]+1][f_idx]^2 - vals_s[f_idx]^2) / vals_s[f_idx]^2 # error measure

                        # calculate the argument of (gsk,hsk)
                        if abs(vals_gs[count_rf[i,1]+1][f_idx]) < eps   #argument is pi/2 or 3pi/2
                            arg = pi/2
                            if abs(vals_hs[count_rf[i,1]+1][f_idx]) < eps
                                @assert vals_s[f_idx] >= 0 "gsk,hsk = 0, S not equal to 0"
                            end
                        elseif abs(vals_hs[count_rf[i,1]+1][f_idx]) < eps   #argument is pi or 0
                            arg = 0
                            if abs(vals_gs[count_rf[i,1]+1][f_idx]) < eps
                                @assert vals_s[f_idx] >= 0 "gsk,hsk = 0, S not equal to 0"
                            end
                        else
                            arg = atan(vals_hs[count_rf[i,1]+1][f_idx]/vals_gs[count_rf[i,1]+1][f_idx])
                        end

                        for tmp_k in count_rf[i,1]+1:K
                            if delta > ref_delta[tmp_k] + eps
                                arg_coef = arg/(pi/2^(tmp_k+1))+1
                                arg_l = Int(round(floor(arg_coef)))     #number of the floor cut in 2^(tmp_k+2) cuts
                                arg_u = Int(round(ceil(arg_coef)))      #number of the ceil cut in 2^(tmp_k+2) cuts

                                idx_l = Int(round((arg_l - 1) * 2^(K-tmp_k)+1))        #number of the floor cut in all cuts
                                idx_u = Int(round((arg_u - 1) * 2^(K-tmp_k)+1))        #number of the ceil cut in all cuts
                                # println("**************point 1**************")
                                if arg_l < arg_u
                                    # @assert outrecord[i,idx_l] == 1 || outrecord[i,idx_u] == 1 "violated with cut added"

                                    if outtensor[i][count_rf[i,1]+1][idx_l] == 1      #the floor cut is not added yet
                                        cut = @build_constraint(
                                            gs[count_rf[i,1]+1][f_idx] * cos((idx_l-1) * pi/2^(K+1)) + hs[count_rf[i,1]+1][f_idx] * sin((idx_l-1) * pi/2^(K+1)) <= s[f_idx]
                                        )
                                        MOI.submit(model, Cons(cb_data), cut)
                                        outtensor[i][count_rf[i,1]+1][idx_l] = 0
                                        updateouttensor(outtensor,num_bra,K)
                                    end

                                    if outtensor[i][count_rf[i,1]+1][idx_u] == 1      #the ceil cut is not added yet
                                        cut = @build_constraint(
                                            gs[count_rf[i,1]+1][f_idx] * cos((idx_u-1) * pi/2^(K+1)) + hs[count_rf[i,1]+1][f_idx] * sin((idx_u-1) * pi/2^(K+1)) <= s[f_idx]
                                        )
                                        MOI.submit(model, Cons(cb_data), cut)
                                        outtensor[i][count_rf[i,1]+1][idx_u] = 0
                                        updateouttensor(outtensor,num_bra,K)
                                    end
                                elseif arg_l == arg_u
                                    # @assert outrecord[i,idx_l] == 1 || outrecord[i,idx_u] == 1 "violated with cut added"
                                    cut = @build_constraint(
                                        gs[count_rf[i,1]+1][f_idx] * cos((idx_l-1) * pi/2^(K+1)) + hs[count_rf[i,1]+1][f_idx] * sin((idx_l-1) * pi/2^(K+1)) <= s[f_idx]
                                    )
                                    MOI.submit(model, Cons(cb_data), cut)
                                    outtensor[i][count_rf[i,1]+1][idx_l] = 0
                                    updateouttensor(outtensor,num_bra,K)
                                end
                                break
                            end
                        end
                    end
                end
            end
        end

        ###--- cone2 ---###
        # inner cut2
        for (i,branch) in ref(pm, n, :branch)
            f_bus = branch["f_bus"]
            t_bus = branch["t_bus"]
            f_idx = (i, f_bus, t_bus)
            tm = branch["tap"]

            # calculate the big M
            bus = ref(pm, n, :bus, f_bus)
            vmi_min = 1/2 * (bus["vmin"]^2 / tm^2 - (branch["rate_a"]^2 * tm^2)/ bus["vmin"]^2)
            vmi_max = bus["vmax"]^2 / (2 *tm^2)
            Mh0 = max(0, -vmi_min)
            hmax = max(vmi_max, -vmi_min)
            M2 = sqrt(hmax^2 + branch["rate_a"]^2)

            for k in K_init+1:K
                # calculate temp_vals[k] based on temp_vals[k-1]
                if k > count_rf[i,2]
                    temp_vals_gi[k+1][f_idx] = temp_vals_gi[k][f_idx] * cos(pi/2^(k+1)) + temp_vals_hi[k][f_idx] * sin(pi/2^(k+1))
                    temp_vals_hi[k+1][f_idx] = abs(-temp_vals_gi[k][f_idx] * sin(pi/2^(k+1)) + temp_vals_hi[k][f_idx] * cos(pi/2^(k+1)))
                end

                if temp_vals_gi[k+1][f_idx]^2 + temp_vals_hi[k+1][f_idx]^2 + eps < vals_vpi[f_idx]^2
                    if !isempty(inrecord_rf[num_bra+i])
                        if inrecord_rf[num_bra+i][end] < K && count_rf[i,2] == K && k == K
                            if vals_vpi[f_idx] * cos(pi/2^(k+2)) > eps + temp_vals_gi[k+1][f_idx] * cos(pi/2^(k+2)) + temp_vals_hi[k+1][f_idx] * sin(pi/2^(k+2))
                                # check if point violate inner cuts 
                                inner_cut = @build_constraint(
                                    vpi[f_idx] * cos(pi/2^(k+2)) <= gi[k+1][f_idx] * cos(pi/2^(k+2)) + hi[k+1][f_idx] * sin(pi/2^(k+2))
                                )
                                MOI.submit(model, Cons(cb_data), inner_cut)
                                append!(inrecord_rf[num_bra+i], k)  # record the history of adding inner cut
                                break   # stop the inner r&f for the current conic surface constraint
                            end
                        end

                        if count_rf[i,2] < k    # check if the k-r&f mapping constraints are added
                            for tmp_k in count_rf[i,2]+1:k
                                # add r&f mapping constraints
                                cons_g = @build_constraint(gi[tmp_k+1][f_idx] == gi[tmp_k][f_idx] * cos(pi/2^(tmp_k+1)) + hi[tmp_k][f_idx] * sin(pi/2^(tmp_k+1)))
                                cons_h1 = @build_constraint(hi_temp[tmp_k][f_idx] == -gi[tmp_k][f_idx] * sin(pi/2^(tmp_k+1)) + hi[tmp_k][f_idx] * cos(pi/2^(tmp_k+1)))
                                cons_h2 = @build_constraint(hi_temp[tmp_k][f_idx] == hi_temp_p[tmp_k][f_idx] - hi_temp_n[tmp_k][f_idx])
                                cons_h3 = @build_constraint(hi[tmp_k+1][f_idx] == hi_temp_p[tmp_k][f_idx] + hi_temp_n[tmp_k][f_idx])
                                cons_h4 = @build_constraint(hi_temp_p[tmp_k][f_idx] <= M2 * flag_hi_temp[tmp_k][f_idx])
                                cons_h5 = @build_constraint(hi_temp_n[tmp_k][f_idx] <= M2 * (1 - flag_hi_temp[tmp_k][f_idx]))
                                cons_bound = @build_constraint(gi[tmp_k+1][f_idx] <= vpi[f_idx])     # this cut is the equivalent cut of the final cuts

                                MOI.submit(model, Cons(cb_data), cons_g)
                                MOI.submit(model, Cons(cb_data), cons_h1)
                                MOI.submit(model, Cons(cb_data), cons_h2)
                                MOI.submit(model, Cons(cb_data), cons_h3)
                                MOI.submit(model, Cons(cb_data), cons_h4)
                                MOI.submit(model, Cons(cb_data), cons_h5)
                                
                                if outer_flag == 5
                                    if outtensor[num_bra+i][tmp_k+1][1] == 1
                                        MOI.submit(model, Cons(cb_data), cons_bound)    #outer cut method 5 needs r&f add the bound outer cut in the same time
                                        outtensor[num_bra+i][tmp_k+1][1] = 0
                                        updateouttensor(outtensor,num_bra,K)
                                    end
                                end

                                count_rf[i,2] = tmp_k   # update the max num of r&f
                            end

                            if inrecord_rf[num_bra+i][end] < k
                                if vals_vpi[f_idx] * cos(pi/2^(k+2)) > eps + temp_vals_gi[k+1][f_idx] * cos(pi/2^(k+2)) + temp_vals_hi[k+1][f_idx] * sin(pi/2^(k+2))
                                    # check if point violate inner cuts 
                                    inner_cut = @build_constraint(
                                        vpi[f_idx] * cos(pi/2^(k+2)) <= gi[k+1][f_idx] * cos(pi/2^(k+2)) + hi[k+1][f_idx] * sin(pi/2^(k+2))
                                    )
                                    MOI.submit(model, Cons(cb_data), inner_cut)
                                    append!(inrecord_rf[num_bra+i], k)  # record the history of adding inner cut
                                    break   # stop the inner r&f for the current conic surface constraint
                                end
                            end
                        end
                    else
                        if k > 0
                            for tmp_k in count_rf[i,2]+1:k
                                # add r&f mapping constraints
                                cons_g = @build_constraint(gi[tmp_k+1][f_idx] == gi[tmp_k][f_idx] * cos(pi/2^(tmp_k+1)) + hi[tmp_k][f_idx] * sin(pi/2^(tmp_k+1)))
                                cons_h1 = @build_constraint(hi_temp[tmp_k][f_idx] == -gi[tmp_k][f_idx] * sin(pi/2^(tmp_k+1)) + hi[tmp_k][f_idx] * cos(pi/2^(tmp_k+1)))
                                cons_h2 = @build_constraint(hi_temp[tmp_k][f_idx] == hi_temp_p[tmp_k][f_idx] - hi_temp_n[tmp_k][f_idx])
                                cons_h3 = @build_constraint(hi[tmp_k+1][f_idx] == hi_temp_p[tmp_k][f_idx] + hi_temp_n[tmp_k][f_idx])
                                cons_h4 = @build_constraint(hi_temp_p[tmp_k][f_idx] <= M2 * flag_hi_temp[tmp_k][f_idx])
                                cons_h5 = @build_constraint(hi_temp_n[tmp_k][f_idx] <= M2 * (1 - flag_hi_temp[tmp_k][f_idx]))
                                cons_bound = @build_constraint(gi[tmp_k+1][f_idx] <= vpi[f_idx])     # this cut is the equivalent cut of the final cuts

                                MOI.submit(model, Cons(cb_data), cons_g)
                                MOI.submit(model, Cons(cb_data), cons_h1)
                                MOI.submit(model, Cons(cb_data), cons_h2)
                                MOI.submit(model, Cons(cb_data), cons_h3)
                                MOI.submit(model, Cons(cb_data), cons_h4)
                                MOI.submit(model, Cons(cb_data), cons_h5)
                               
                                if outer_flag == 5
                                    if outtensor[num_bra+i][tmp_k+1][1] == 1
                                        MOI.submit(model, Cons(cb_data), cons_bound)    #outer cut method 5 needs r&f add the bound outer cut in the same time
                                        outtensor[num_bra+i][tmp_k+1][1] = 0
                                        updateouttensor(outtensor,num_bra,K)
                                    end
                                end

                                count_rf[i,2] = tmp_k   # update the max num of r&f
                            end
                        end
                        # check if 0-inner cuts are violated? if so, add the inner cuts
                        if vals_vpi[f_idx] * cos(pi/2^(k+2)) > eps + temp_vals_gi[k+1][f_idx] * cos(pi/2^(k+2)) + temp_vals_hi[k+1][f_idx] * sin(pi/2^(k+2))
                            # check if point violate inner cuts 
                            inner_cut = @build_constraint(
                                vpi[f_idx] * cos(pi/2^(k+2)) <= gi[k+1][f_idx] * cos(pi/2^(k+2)) + hi[k+1][f_idx] * sin(pi/2^(k+2))
                            )
                            MOI.submit(model, Cons(cb_data), inner_cut)
                            append!(inrecord_rf[num_bra+i], k)  # record the history of adding inner cut
                            break   # stop the inner r&f for the current conic surface constraint
                        end
                    end
                end
            end
        end

        # outer cut2
        if constraints_flag == 4
            if outer_flag == 1
                # outer cut2 ::method 1: add outer cut in r&f way
                for (i,branch) in ref(pm, n, :branch)
                    f_bus = branch["f_bus"]
                    t_bus = branch["t_bus"]
                    f_idx = (i, f_bus, t_bus)
                    tm = branch["tap"]
        
                    # calculate the big M
                    bus = ref(pm, n, :bus, f_bus)
                    vmi_min = 1/2 * (bus["vmin"]^2 / tm^2 - (branch["rate_a"]^2 * tm^2)/ bus["vmin"]^2)
                    vmi_max = bus["vmax"]^2 / (2 *tm^2)
                    Mh0 = max(0, -vmi_min)
                    hmax = max(vmi_max, -vmi_min)
                    M2 = sqrt(hmax^2 + branch["rate_a"]^2)

                    for k in K_init:K
                        if k > count_rf[i,2]
                            # calculate temp_vals[k] based on temp_vals[k-1]
                            temp_vals_gi[k+1][f_idx] = temp_vals_gi[k][f_idx] * cos(pi/2^(k+1)) + temp_vals_hi[k][f_idx] * sin(pi/2^(k+1))
                            temp_vals_hi[k+1][f_idx] = abs(-temp_vals_gi[k][f_idx] * sin(pi/2^(k+1)) + temp_vals_hi[k][f_idx] * cos(pi/2^(k+1)))
                        end

                        if k == count_rf[i,2] # need not to add r&f mapping and may need add outer cut 
                            if k == K_init || k == outrecord_rf[num_bra+i][end]
                                continue        # k-outer cut already added
                            elseif (temp_vals_gi[k+1][f_idx] > vals_vpi[f_idx] + eps) || 
                                (temp_vals_gi[k+1][f_idx] * cos(pi/2^(k+1)) + temp_vals_hi[k+1][f_idx] * sin(pi/2^(k+1)) > vals_vpi[f_idx] + eps)
                                # k-outer cut is violated
                                outer_cut1 = @build_constraint(gi[k+1][f_idx] <= vpi[f_idx])
                                outer_cut2 = @build_constraint(gi[k+1][f_idx] * cos(pi/2^(k+1)) + hi[k+1][f_idx] * sin(pi/2^(k+1)) <= vpi[f_idx])

                                MOI.submit(model, Cons(cb_data), outer_cut1)
                                MOI.submit(model, Cons(cb_data), outer_cut2)

                                append!(outrecord_rf[num_bra+i], k)  # record the history of adding outer cut
                                break   # stop the outer r&f for the current conic surface constraint
                            else
                                continue       # k-outer cut already added
                            end

                        end

                        if temp_vals_gi[k+1][f_idx]^2 + temp_vals_hi[k+1][f_idx]^2 > vals_vpi[f_idx]^2 + eps  # check if point is located outside the cone 
                            if outrecord_rf[num_bra+i][end] < k && count_rf[i,2] == K && k == K      # all r&f mappings are added but the outer cut is not added
                                if (temp_vals_gi[k+1][f_idx] > vals_vpi[f_idx] + eps) || 
                                    (temp_vals_gi[k+1][f_idx] * cos(pi/2^(k+1)) + temp_vals_hi[k+1][f_idx] * sin(pi/2^(k+1)) > vals_vpi[f_idx] + eps)
                                    # k-outer cut is violated
                                    outer_cut1 = @build_constraint(gi[k+1][f_idx] <= vpi[f_idx])
                                    outer_cut2 = @build_constraint(gi[k+1][f_idx] * cos(pi/2^(k+1)) + hi[k+1][f_idx] * sin(pi/2^(k+1)) <= vpi[f_idx])

                                    MOI.submit(model, Cons(cb_data), outer_cut1)
                                    MOI.submit(model, Cons(cb_data), outer_cut2)

                                    append!(outrecord_rf[num_bra+i], k)  # record the history of adding outer cut
                                    break   # stop the outer r&f for the current conic surface constraint
                                end
                            end

                            if count_rf[i,2] < k    # check if the k-r&f mapping constraints are added
                                for tmp_k in count_rf[i,2]+1:k
                                    # add r&f mapping constraints
                                    cons_g = @build_constraint(gi[tmp_k+1][f_idx] == gi[tmp_k][f_idx] * cos(pi/2^(tmp_k+1)) + hi[tmp_k][f_idx] * sin(pi/2^(tmp_k+1)))
                                    cons_h1 = @build_constraint(hi_temp[tmp_k][f_idx] == -gi[tmp_k][f_idx] * sin(pi/2^(tmp_k+1)) + hi[tmp_k][f_idx] * cos(pi/2^(tmp_k+1)))
                                    cons_h2 = @build_constraint(hi_temp[tmp_k][f_idx] == hi_temp_p[tmp_k][f_idx] - hi_temp_n[tmp_k][f_idx])
                                    cons_h3 = @build_constraint(hi[tmp_k+1][f_idx] == hi_temp_p[tmp_k][f_idx] + hi_temp_n[tmp_k][f_idx])
                                    cons_h4 = @build_constraint(hi_temp_p[tmp_k][f_idx] <= M2 * flag_hi_temp[tmp_k][f_idx])
                                    cons_h5 = @build_constraint(hi_temp_n[tmp_k][f_idx] <= M2 * (1 - flag_hi_temp[tmp_k][f_idx]))

                                    MOI.submit(model, Cons(cb_data), cons_g)
                                    MOI.submit(model, Cons(cb_data), cons_h1)
                                    MOI.submit(model, Cons(cb_data), cons_h2)
                                    MOI.submit(model, Cons(cb_data), cons_h3)
                                    MOI.submit(model, Cons(cb_data), cons_h4)
                                    MOI.submit(model, Cons(cb_data), cons_h5)
                                
                                    count_rf[i,2] = tmp_k   # update the max num of r&f
                                end

                                # check if temp_vals[k+1] violate the inner cut. if violated, add the cut and break
                                if (temp_vals_gi[k+1][f_idx] > vals_vpi[f_idx] + eps) || 
                                    (temp_vals_gi[k+1][f_idx] * cos(pi/2^(k+1)) + temp_vals_hi[k+1][f_idx] * sin(pi/2^(k+1)) > vals_vpi[f_idx] + eps)
                                    # k-outer cut is violated
                                    outer_cut1 = @build_constraint(gi[k+1][f_idx] <= vpi[f_idx])
                                    outer_cut2 = @build_constraint(gi[k+1][f_idx] * cos(pi/2^(k+1)) + hi[k+1][f_idx] * sin(pi/2^(k+1)) <= vpi[f_idx])

                                    MOI.submit(model, Cons(cb_data), outer_cut1)
                                    MOI.submit(model, Cons(cb_data), outer_cut2)

                                    append!(outrecord_rf[num_bra+i], k)  # record the history of adding outer cut
                                    break   # stop the outer r&f for the current conic surface constraint
                                end
                            end
                        end

                    end
                end
            elseif outer_flag == 2
                # outer cut2 ::method 2: add central cut when violation occurs
                for (i,branch) in ref(pm, n, :branch)
                    f_bus = branch["f_bus"]
                    t_bus = branch["t_bus"]
                    f_idx = (i, f_bus, t_bus)
                    tm = branch["tap"]

                    # calculate the big M
                    bus = ref(pm, n, :bus, f_bus)
                    vmi_min = 1/2 * (bus["vmin"]^2 / tm^2 - (branch["rate_a"]^2 * tm^2)/ bus["vmin"]^2)
                    vmi_max = bus["vmax"]^2 / (2 *tm^2)
                    Mh0 = max(0, -vmi_min)
                    hmax = max(vmi_max, -vmi_min)
                    M2 = sqrt(hmax^2 + branch["rate_a"]^2)

                    for k in count_rf[i,2]:K
                        if outrecord_rf[i+num_bra][end] >= K    # number of outer cuts reached limit
                            break
                        end

                        # calculate temp_vals[k] based on temp_vals[k-1]
                        if k > count_rf[i,2]
                            temp_vals_gi[k+1][f_idx] = temp_vals_gi[k][f_idx] * cos(pi/2^(k+1)) + temp_vals_hi[k][f_idx] * sin(pi/2^(k+1))
                            temp_vals_hi[k+1][f_idx] = abs(-temp_vals_gi[k][f_idx] * sin(pi/2^(k+1)) + temp_vals_hi[k][f_idx] * cos(pi/2^(k+1)))
                        end

                        if temp_vals_gi[k+1][f_idx]^2 + temp_vals_hi[k+1][f_idx]^2 > vals_vpi[f_idx]^2 + eps  # check if point is located outside the cone 
                            if k > count_rf[i,2]
                                cons_g = @build_constraint(gi[k+1][f_idx] == gi[k][f_idx] * cos(pi/2^(k+1)) + hi[k][f_idx] * sin(pi/2^(k+1)))
                                cons_h1 = @build_constraint(hi_temp[k][f_idx] == -gi[k][f_idx] * sin(pi/2^(k+1)) + hi[k][f_idx] * cos(pi/2^(k+1)))
                                cons_h2 = @build_constraint(hi_temp[k][f_idx] == hi_temp_p[k][f_idx] - hi_temp_n[k][f_idx])
                                cons_h3 = @build_constraint(hi[k+1][f_idx] == hi_temp_p[k][f_idx] + hi_temp_n[k][f_idx])
                                cons_h4 = @build_constraint(hi_temp_p[k][f_idx] <= M2 * flag_hi_temp[k][f_idx])
                                cons_h5 = @build_constraint(hi_temp_n[k][f_idx] <= M2 * (1 - flag_hi_temp[k][f_idx]))
                                cons_bound = @build_constraint(gi[k+1][f_idx] <= vpi[f_idx])     # this cut is the equivalent cut of the final cuts

                                MOI.submit(model, Cons(cb_data), cons_g)
                                MOI.submit(model, Cons(cb_data), cons_h1)
                                MOI.submit(model, Cons(cb_data), cons_h2)
                                MOI.submit(model, Cons(cb_data), cons_h3)
                                MOI.submit(model, Cons(cb_data), cons_h4)
                                MOI.submit(model, Cons(cb_data), cons_h5)
                                MOI.submit(model, Cons(cb_data), cons_bound)    #outer cut method 2 needs r&f add the bound outer cut in the same time

                                count_rf[i,2] = k   # update the max num of r&f
                                append!(outrecord_rf[num_bra+i], k)  # record the history of adding outer cut

                            end

                            if k == count_rf[i,2]
                                tmp_delta = (temp_vals_gi[k+1][f_idx]^2 + temp_vals_hi[k+1][f_idx]^2 - vals_vpi[f_idx]^2) / vals_vpi[f_idx]^2 # error measure

                                # prepare the outer cuts
                                cons_tan1 = @build_constraint(gi[k+1][f_idx] * cos(pi/2^(k+2)) + hi[k+1][f_idx] * sin(pi/2^(k+2)) <= vpi[f_idx])
                                cons_tan2 = @build_constraint(gi[k+1][f_idx] * cos(pi/2^(k+3)) + hi[k+1][f_idx] * sin(pi/2^(k+3)) <= vpi[f_idx])
                                cons_tan3 = @build_constraint(gi[k+1][f_idx] * cos(3pi/2^(k+3)) + hi[k+1][f_idx] * sin(3pi/2^(k+3)) <= vpi[f_idx])
                                cons_tan4 = @build_constraint(gi[k+1][f_idx] * cos(pi/2^(k+4)) + hi[k+1][f_idx] * sin(pi/2^(k+4)) <= vpi[f_idx])
                                cons_tan5 = @build_constraint(gi[k+1][f_idx] * cos(3pi/2^(k+4)) + hi[k+1][f_idx] * sin(3pi/2^(k+4)) <= vpi[f_idx])
                                cons_tan6 = @build_constraint(gi[k+1][f_idx] * cos(5pi/2^(k+4)) + hi[k+1][f_idx] * sin(5pi/2^(k+4)) <= vpi[f_idx])
                                cons_tan7 = @build_constraint(gi[k+1][f_idx] * cos(7pi/2^(k+4)) + hi[k+1][f_idx] * sin(7pi/2^(k+4)) <= vpi[f_idx])

                                if tmp_delta > ref_delta[k+1]+1e-10
                                     # central cut should be add
                                    # @assert outrecord_rf[i+num_bra][end]<k+1 "Caution: out1"       # the central cut should not be added yet

                                    if k+1 > K     
                                        break
                                    else
                                        MOI.submit(model, Cons(cb_data), cons_tan1) 
                                        append!(outrecord_rf[i+num_bra], k+1)        # record the history of adding outer cut 
                                        break
                                    end
                                elseif tmp_delta > ref_delta[k+2]+1e-10   # 2 subcentral cuts should be add
                                    # @assert outrecord_rf[i+num_bra][end]<k+2 "Caution: out2"       # the subcentral cuts should not be added yet

                                    if k+2 > K     
                                        break
                                    else
                                        MOI.submit(model, Cons(cb_data), cons_tan1) 
                                        MOI.submit(model, Cons(cb_data), cons_tan2) 
                                        MOI.submit(model, Cons(cb_data), cons_tan3) 
                                        append!(outrecord_rf[i+num_bra], k+2)        # record the history of adding outer cut 
                                        break
                                    end
                                elseif tmp_delta > ref_delta[k+3]+1e-10   # 4 subcentral cuts should be add
                                    # @assert outrecord_rf[i+num_bra][end]<k+3 "Caution: out3"       # the subcentral cuts should not be added yet

                                    if k+3 > K     
                                        break
                                    else
                                        MOI.submit(model, Cons(cb_data), cons_tan1) 
                                        MOI.submit(model, Cons(cb_data), cons_tan2) 
                                        MOI.submit(model, Cons(cb_data), cons_tan3) 
                                        MOI.submit(model, Cons(cb_data), cons_tan4) 
                                        MOI.submit(model, Cons(cb_data), cons_tan5) 
                                        MOI.submit(model, Cons(cb_data), cons_tan6)
                                        MOI.submit(model, Cons(cb_data), cons_tan7)
                                        append!(outrecord_rf[i+num_bra], k+3)        # record the history of adding inner cut 
                                        break
                                    end
                                else 
                                    if k >= K-3
                                        break   
                                    else 
                                        continue
                                    end
                                end
                            end
                        end
                    end
                end
            elseif outer_flag == 3
                # # outer cut2 ::method 3: add cuts based on arguments of (Sr/baseMVA,VrI_add)
                # calculate delta
                ref_delta = zeros(K+3)
                for k in 1:K+3
                    ref_delta[k]=(1/(cos(pi/2^(k+2))).^2)-1  
                end

                vals_s = callback_value.(cb_data, s)
                vals_vpi = callback_value.(cb_data, vpi)
                vals_vmi = callback_value.(cb_data, vmi)

                for (i,branch) in ref(pm, n, :branch)
                    f_bus = branch["f_bus"]
                    t_bus = branch["t_bus"]
                    f_idx = (i, f_bus, t_bus)
                    tm = branch["tap"]

                    # calculate the big M
                    bus = ref(pm, n, :bus, f_bus)
                    vmi_min = 1/2 * (bus["vmin"]^2 / tm^2 - (branch["rate_a"]^2 * tm^2)/ bus["vmin"]^2)
                    vmi_max = bus["vmax"]^2 / (2 *tm^2)
                    Mh0 = max(0, -vmi_min)
                    hmax = max(vmi_max, -vmi_min)
                    M2 = sqrt(hmax^2 + branch["rate_a"]^2)

                    if vals_s[f_idx]^2 + vals_vmi[f_idx]^2 > eps + vals_vpi[f_idx]^2
                        delta = (vals_s[f_idx]^2 + vals_vmi[f_idx]^2 - vals_vpi[f_idx]^2) / vals_vpi[f_idx]^2 # error measure

                        # calculate the argument of (gsk,hsk)
                        if vals_s[f_idx] == 0
                            if vals_vmi[f_idx] > 0
                                arg = pi/2
                            elseif vals_vmi[f_idx] < 0
                                arg = 3pi/2
                            else
                                @assert vals_vpi[f_idx] != 0 "gsk,hsk = 0, S not equal to 0"
                            end
                        elseif vals_vmi[f_idx] == 0
                            if vals_s[f_idx] > 0
                                arg = 0
                            elseif vals_s[f_idx] < 0
                                arg = pi
                            else
                                @assert vals_vpi[f_idx] != 0 "gsk,hsk = 0, S not equal to 0"
                            end
                        else
                            if vals_s[f_idx] > 0 && vals_vmi[f_idx] > 0
                                arg = atan(vals_vmi[f_idx]/vals_s[f_idx])
                            elseif vals_s[f_idx] < 0 && vals_vmi[f_idx] > 0
                                arg = pi + atan(vals_vmi[f_idx]/vals_s[f_idx])
                            elseif vals_s[f_idx] < 0 && vals_vmi[f_idx] < 0
                                arg = pi + atan(vals_vmi[f_idx]/vals_s[f_idx])
                            elseif vals_s[f_idx] > 0 && vals_vmi[f_idx] < 0
                                arg = 2pi + atan(vals_vmi[f_idx]/vals_s[f_idx])
                            end
                        end

                        for tmp_k in 1:K
                            if delta > ref_delta[tmp_k] + eps
                                arg_coef = arg/(pi/2^(tmp_k+1))+1
                                arg_l = Int(round(floor(arg_coef)))     #number of the floor cut in 2^(tmp_k+2) cuts
                                arg_u = Int(round(ceil(arg_coef)))      #number of the ceil cut in 2^(tmp_k+2) cuts

                                if arg_u == 2^(tmp_k+2)+1
                                    if arg_l < arg_u
                                        arg_u = 1
                                    elseif arg_l == arg_u
                                        arg_u = 1
                                        arg_l = 1
                                    end
                                end

                                idx_l = Int(round((arg_l - 1) * 2^(K-tmp_k)+1))        #number of the floor cut in all cuts
                                idx_u = Int(round((arg_u - 1) * 2^(K-tmp_k)+1))        #number of the ceil cut in all cuts

                                if arg_l < arg_u || (arg_u == 1 && arg_l != 1)                 #solution is not located at the border

                                    # @assert outrecord_pq[i+num_bra,idx_l] == 1 || outrecord_pq[i+num_bra,idx_u] == 1 "violated with cut added"

                                    if outrecord_pq[i+num_bra,idx_l] == 1      #the floor cut is not added yet
                                        cut = @build_constraint(
                                            s[f_idx] * cos((idx_l-1) * pi/2^(K+1)) + vmi[f_idx] * sin((idx_l-1) * pi/2^(K+1)) <= vpi[f_idx]
                                        )
                                        MOI.submit(model, Cons(cb_data), cut)
                                        outrecord_pq[i+num_bra,idx_l] = 0
                                    end
                                    if outrecord_pq[i+num_bra,idx_u] == 1      #the ceil cut is not added yet
                                        cut = @build_constraint(
                                            s[f_idx] * cos((idx_u-1) * pi/2^(K+1)) + vmi[f_idx] * sin((idx_u-1) * pi/2^(K+1)) <= vpi[f_idx]
                                        )
                                        MOI.submit(model, Cons(cb_data), cut)
                                        outrecord_pq[i+num_bra,idx_u] = 0
                                    end
                                elseif arg_l == arg_u || (arg_u == 1 && arg_l == 1)
                                    # @assert outrecord_pq[i+num_bra,idx_l] == 1 || outrecord_pq[i+num_bra,idx_u] == 1 "violated with cut added"
                                    if outrecord_pq[i+num_bra,idx_l] == 1 
                                        cut = @build_constraint(
                                            s[f_idx] * cos((idx_l-1) * pi/2^(K+1)) + vmi[f_idx] * sin((idx_l-1) * pi/2^(K+1)) <= vpi[f_idx]
                                        )
                                        MOI.submit(model, Cons(cb_data), cut)
                                        outrecord_pq[i+num_bra,idx_l] = 0
                                    end
                                end
                                break
                            end
                        end
                    end
                end
            elseif outer_flag == 4
                # outer cut2 ::method 4: add cuts based on arguments of (gi0,hi0)
                # calculate delta
                ref_delta = zeros(Float64,K)
                for k in 1:K
                    ref_delta[k]=(1/(cos(pi/2^(k+2))).^2)-1
                end

                for (i,branch) in ref(pm, n, :branch)
                    f_bus = branch["f_bus"]
                    t_bus = branch["t_bus"]
                    f_idx = (i, f_bus, t_bus)
                    tm = branch["tap"]

                    # calculate the big M
                    bus = ref(pm, n, :bus, f_bus)
                    vmi_min = 1/2 * (bus["vmin"]^2 / tm^2 - (branch["rate_a"]^2 * tm^2)/ bus["vmin"]^2)
                    vmi_max = bus["vmax"]^2 / (2 *tm^2)
                    Mh0 = max(0, -vmi_min)
                    hmax = max(vmi_max, -vmi_min)
                    M2 = sqrt(hmax^2 + branch["rate_a"]^2)

                    if vals_gi[1][f_idx]^2 + vals_hi[1][f_idx]^2 > eps + vals_vpi[f_idx]^2
                        delta = (vals_gi[1][f_idx]^2 + vals_hi[1][f_idx]^2 - vals_vpi[f_idx]^2) / vals_vpi[f_idx]^2 # error measure

                        # calculate the argument of (gsk,hsk)
                        if vals_gi[1][f_idx] == 0
                            arg = pi/2
                            if vals_hi[1][f_idx] == 0
                                @assert vals_vpi[f_idx] != 0 "gs1,hs1 = 0, S not equal to 0"
                            end
                        elseif vals_hi[1][f_idx] == 0
                            arg = 0
                            if vals_gi[1][f_idx] == 0
                                @assert vals_vpi[f_idx] != 0 "gs1,hs1 = 0, S not equal to 0"
                            end
                        else
                            arg = atan(vals_hi[1][f_idx]/vals_gi[1][f_idx])
                        end

                        for tmp_k in 1:K
                            if delta > ref_delta[tmp_k] + eps
                                arg_coef = arg/(pi/2^(tmp_k+1))+1
                                arg_l = Int(round(floor(arg_coef)))        #number of the floor cut in 2^(tmp_k+2) cuts
                                arg_u = Int(round(ceil(arg_coef)))         #number of the ceil cut in 2^(tmp_k+2) cuts
        
                                idx_l = Int(round((arg_l - 1) * 2^(K-tmp_k)+1))        #number of the floor cut in all cuts
                                idx_u = Int(round((arg_u - 1) * 2^(K-tmp_k)+1))        #number of the ceil cut in all cuts
        
                                if arg_l < arg_u
                                    # @assert outrecord[i+num_bra,idx_l] == 1 || outrecord[i+num_bra,idx_u] == 1 "violated with cut added"
        
                                    if outrecord[i+num_bra,idx_l] == 1      #the floor cut is not added yet
                                        cut = @build_constraint(
                                            gi[1][f_idx] * cos((idx_l-1) * pi/2^(K+1)) + hi[1][f_idx] * sin((idx_l-1) * pi/2^(K+1)) <= vpi[f_idx]
                                        )
                                        MOI.submit(model, Cons(cb_data), cut)
                                        outrecord[i+num_bra,idx_l] = 0
                                    end

                                    if outrecord[i+num_bra,idx_u] == 1      #the ceil cut is not added yet
                                        cut = @build_constraint(
                                            gi[1][f_idx] * cos((idx_u-1) * pi/2^(K+1)) + hi[1][f_idx] * sin((idx_u-1) * pi/2^(K+1)) <= vpi[f_idx]
                                        )
                                        MOI.submit(model, Cons(cb_data), cut)
                                        outrecord[i+num_bra,idx_u] = 0
                                    end
                                elseif arg_l == arg_u
                                    # @assert outrecord[i+num_bra,idx_l] == 1 || outrecord[i+num_bra,idx_u] == 1 "violated with cut added"
                                    cut = @build_constraint(
                                        gi[1][f_idx] * cos((idx_l-1) * pi/2^(K+1)) + hi[1][f_idx] * sin((idx_l-1) * pi/2^(K+1)) <= vpi[f_idx]
                                    )
                                    MOI.submit(model, Cons(cb_data), cut)
                                    outrecord[i+num_bra,idx_l] = 0
                                end
                                break
                            end
                        end
                    end
                end
            elseif outer_flag == 5
                # outer cut1 ::method 5: add cuts based on arguments of (gik,hik) where k = count_rf
                # calculate delta
                ref_delta = zeros(Float64,K)
                for k in 1:K
                    ref_delta[k]=(1/(cos(pi/2^(k+2))).^2)-1
                end

                for (i,branch) in ref(pm, n, :branch)
                    f_bus = branch["f_bus"]
                    t_bus = branch["t_bus"]
                    f_idx = (i, f_bus, t_bus)
                    tm = branch["tap"]

                    # calculate the big M
                    bus = ref(pm, n, :bus, f_bus)
                    vmi_min = 1/2 * (bus["vmin"]^2 / tm^2 - (branch["rate_a"]^2 * tm^2)/ bus["vmin"]^2)
                    vmi_max = bus["vmax"]^2 / (2 *tm^2)
                    Mh0 = max(0, -vmi_min)
                    hmax = max(vmi_max, -vmi_min)
                    M2 = sqrt(hmax^2 + branch["rate_a"]^2)

                    if vals_gi[count_rf[i,2]+1][f_idx]^2 + vals_hi[count_rf[i,2]+1][f_idx]^2 > eps + vals_vpi[f_idx]^2
                        delta = (vals_gi[count_rf[i,2]+1][f_idx]^2 + vals_hi[count_rf[i,2]+1][f_idx]^2 - vals_vpi[f_idx]^2) / vals_vpi[f_idx]^2 # error measure

                        # calculate the argument of (gsk,hsk)
                        if abs(vals_gi[count_rf[i,2]+1][f_idx]) < eps
                            arg = pi/2
                            if abs(vals_hi[count_rf[i,2]+1][f_idx]) < eps
                                @assert abs(vals_vpi[f_idx]) >= 0 "gsk,hsk = 0, S not equal to 0"
                            end
                        elseif abs(vals_hi[count_rf[i,2]+1][f_idx]) < eps
                            arg = 0
                            if abs(vals_gi[count_rf[i,2]+1][f_idx]) < eps
                                @assert abs(vals_vpi[f_idx]) >= 0 "gsk,hsk = 0, S not equal to 0"
                            end
                        else
                            arg = atan(vals_hi[count_rf[i,2]+1][f_idx]/vals_gi[count_rf[i,2]+1][f_idx])
                        end
                        for tmp_k in count_rf[i,2]+1:K
                            if delta > ref_delta[tmp_k] + eps
                                arg_coef = arg/(pi/2^(tmp_k+1))+1
                                arg_l = Int(round(floor(arg_coef)))     #number of the floor cut in 2^(tmp_k+2) cuts
                                arg_u = Int(round(ceil(arg_coef)))      #number of the ceil cut in 2^(tmp_k+2) cuts

                                idx_l = Int(round((arg_l - 1) * 2^(K-tmp_k)+1))        #number of the floor cut in all cuts
                                idx_u = Int(round((arg_u - 1) * 2^(K-tmp_k)+1))        #number of the ceil cut in all cuts

                                if arg_l < arg_u
                                    # @assert outrecord[i,idx_l] == 1 || outrecord[i,idx_u] == 1 "violated with cut added"

                                    if outtensor[num_bra+i][count_rf[i,2]+1][idx_l] == 1      #the floor cut is not added yet
                                        cut = @build_constraint(
                                            gi[count_rf[i,2]+1][f_idx] * cos((idx_l-1) * pi/2^(K+1)) + hi[count_rf[i,2]+1][f_idx] * sin((idx_l-1) * pi/2^(K+1)) <= vpi[f_idx]
                                        )
                                        MOI.submit(model, Cons(cb_data), cut)
                                        outtensor[num_bra+i][count_rf[i,2]+1][idx_l] = 0
                                        updateouttensor(outtensor,num_bra,K)
                                    end

                                    if outtensor[num_bra+i][count_rf[i,2]+1][idx_u] == 1      #the ceil cut is not added yet
                                        cut = @build_constraint(
                                            gi[count_rf[i,2]+1][f_idx] * cos((idx_u-1) * pi/2^(K+1)) + hi[count_rf[i,2]+1][f_idx] * sin((idx_u-1) * pi/2^(K+1)) <= vpi[f_idx]
                                        )
                                        MOI.submit(model, Cons(cb_data), cut)
                                        outtensor[num_bra+i][count_rf[i,2]+1][idx_u] = 0
                                        updateouttensor(outtensor,num_bra,K)
                                    end
                                elseif arg_l == arg_u
                                    # @assert outrecord[i,idx_l] == 1 || outrecord[i,idx_u] == 1 "violated with cut added"
                                    cut = @build_constraint(
                                        gi[count_rf[i,2]+1][f_idx] * cos((idx_l-1) * pi/2^(K+1)) + hi[count_rf[i,2]+1][f_idx] * sin((idx_l-1) * pi/2^(K+1)) <= vpi[f_idx]
                                    )
                                    MOI.submit(model, Cons(cb_data), cut)
                                    outtensor[num_bra+i][count_rf[i,2]+1][idx_l] = 0
                                    updateouttensor(outtensor,num_bra,K)
                                end
                                break
                            end
                        end
                    end
                end
            end
        end
    end
end

function assign_warm_start(pm::AbstractSOCDRBFModel, warm_res, warm_start_model, n=nw_id_default)
    model = pm.model

    DEBUG = false
    # DEBUG = true
    ccm = var(pm, n, :ccm)
    p = var(pm, n, :p)
    q = var(pm, n, :q)
    pg = var(pm, n, :pg)
    qg = var(pm, n, :qg)
    w = var(pm, n, :w)
    
    sol = warm_res["solution"]
    if warm_start_model === IVRPowerModel
        # set values for branch variables
        for (i,branch) in ref(pm, n, :branch)
            f_bus = branch["f_bus"]
            t_bus = branch["t_bus"]
            f_idx = (i, f_bus, t_bus)
            t_idx = (i, t_bus, f_bus)
            tm = branch["tap"]

            # set values for ccm
            val_ccm = sol["branch"][string(i)]["csr_fr"]^2 + sol["branch"][string(i)]["csi_fr"]^2
            if i == 1
                println("branch $i: ccm = $val_ccm")
            end
            JuMP.set_start_value(ccm[i], val_ccm)
            if DEBUG
               JuMP.fix(ccm[i], val_ccm, force=true) 
            end

            val_pf = sol["branch"][string(i)]["pf"]
            val_qf = sol["branch"][string(i)]["qf"]
            if i == 1
                println("branch $i: pf = $val_pf, qf = $val_qf")
            end
            JuMP.set_start_value(p[f_idx], val_pf)
            JuMP.set_start_value(q[f_idx], val_qf)
            if DEBUG
                JuMP.fix(p[f_idx], val_pf, force=true)
                JuMP.fix(q[f_idx], val_qf, force=true)
            end

            # val_pt = sol["branch"][string(i)]["pt"]
            # val_qt = sol["branch"][string(i)]["qt"]
            # if i == 1
            #     println("branch $i: pt = $val_pt, qt = $val_qt")
            # end
            # JuMP.set_start_value(p[t_idx], val_pt)
            # JuMP.set_start_value(q[t_idx], val_qt)
            # JuMP.fix(p[t_idx], val_pt, force=true)
            # JuMP.fix(q[t_idx], val_qt, force=true)
        end
        # set values for gen variables
        for (i,gen) in ref(pm, n, :gen)
            val_pg = sol["gen"][string(i)]["pg"]
            val_qg = sol["gen"][string(i)]["qg"]
            JuMP.set_start_value(pg[i], val_pg)
            JuMP.set_start_value(qg[i], val_qg)
            if DEBUG
                JuMP.fix(pg[i], val_pg, force=true)
                JuMP.fix(qg[i], val_qg, force=true)
            end
        end
        # set values for bus variables
        for (i,bus) in ref(pm, n, :bus)
            val_vr = sol["bus"][string(i)]["vr"]
            val_vi = sol["bus"][string(i)]["vi"]
            val_w = val_vr^2 + val_vi^2
            JuMP.set_start_value(w[i], val_w)
            if DEBUG
                JuMP.fix(w[i], val_w, force=true)
            end
        end
    elseif warm_start_model === ACPPowerModel
        # set values for branch variables
        for (i,branch) in ref(pm, n, :branch)
            f_bus = branch["f_bus"]
            t_bus = branch["t_bus"]
            f_idx = (i, f_bus, t_bus)
            t_idx = (i, t_bus, f_bus)
            tm = branch["tap"]

            # set values for power flow variables
            val_pf = sol["branch"][string(i)]["pf"]
            val_qf = sol["branch"][string(i)]["qf"]
            if i == 1
                println("branch $i: pf = $val_pf, qf = $val_qf")
            end
            JuMP.set_start_value(p[f_idx], val_pf)
            JuMP.set_start_value(q[f_idx], val_qf)
            if DEBUG
                JuMP.fix(p[f_idx], val_pf, force=true)
                JuMP.fix(q[f_idx], val_qf, force=true)
            end

            # calculate values of ccm
            val_ccm = tm^2 * (val_pf^2 + val_qf^2) / sol["bus"][string(f_bus)]["vm"]^2
            if i == 1
                println("branch $i: ccm = $val_ccm")
            end
            JuMP.set_start_value(ccm[i], val_ccm)
            if DEBUG
                JuMP.fix(ccm[i], val_ccm, force=true)
            end

        end
        # set values for gen variables
        for (i,gen) in ref(pm, n, :gen)
            val_pg = sol["gen"][string(i)]["pg"]
            val_qg = sol["gen"][string(i)]["qg"]
            JuMP.set_start_value(pg[i], val_pg)
            JuMP.set_start_value(qg[i], val_qg)
            if DEBUG
                JuMP.fix(pg[i], val_pg, force=true)
                JuMP.fix(qg[i], val_qg, force=true)
            end
        end
        # set values for bus variables
        for (i,bus) in ref(pm, n, :bus)
            val_vm = sol["bus"][string(i)]["vm"]
            
            val_w = val_vm^2
            JuMP.set_start_value(w[i], val_w)
            if DEBUG
                JuMP.fix(w[i], val_w, force=true)
            end
        end
    else
        #TODO: add warm start for other models
    end

    
end