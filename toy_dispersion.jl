using Plots

# An analytical dispersion relation
function dispersion_relation(J, S, H, K, L, Δ)
    return J * S * (3 - cos(2π * H) - cos(2π * K) - cos(2π * L)) + Δ
end

# Create qpts, a 3D array where each element contains (H, K, L)
function generate_qpts(H_vals, K_vals, L_vals)
    qpts = Array{Tuple{Float64, Float64, Float64}}(undef, length(H_vals), length(K_vals), length(L_vals))
    
    for i in 1:length(H_vals)
        for j in 1:length(K_vals)
            for k in 1:length(L_vals)
                qpts[i, j, k] = (H_vals[i], K_vals[j], L_vals[k])
            end
        end
    end
    return qpts
end

# Calculate the dispersion modes using the qpts object and store both modes in a single array
function calculate_modes_from_qpts(qpts)
    # We will store 3 modes for demonstration (can be any number of modes)
    modes = Array{NTuple{3, Float64}}(undef, size(qpts)...)  # Each element will store multiple modes
    
    for i in 1:size(qpts, 1)
        for j in 1:size(qpts, 2)
            for k in 1:size(qpts, 3)
                H, K, L = qpts[i, j, k]
                mode1 = dispersion_relation(5, 1, H, K, L, 0)  # Mode 1: S = 1, J = 5, Δ = 0
                mode2 = dispersion_relation(2, 1, H, K, L, 5)  # Mode 2: S = 1, J = 2, Δ = 5
                mode3 = dispersion_relation(3, 1, H, K, L, 2)  # Mode 3: S = 1, J = 3, Δ = 2
                modes[i, j, k] = (mode1, mode2, mode3)  # Store all modes as a tuple
            end
        end
    end
    return modes
end

# Perform 1D slicing delta-function along H, K, or L
function oneD_mode_slice(modes, qpts, slice_definition)
    h_def, k_def, l_def = slice_definition
    
    # Determine which dimension is a range and which are fixed
    if typeof(h_def) <: Tuple
        # We're slicing along H (h_def is a range)
        H_vals = [qpts[i, 1, 1][1] for i in 1:size(qpts, 1)]  # Extract H values from qpts
        H_range = filter(h -> h >= h_def[1] && h <= h_def[2], H_vals)  # Create H range within the slice range
        K_idx = findfirst(x -> qpts[1, x, 1][2] >= k_def, 1:size(qpts, 2))  # Find K index
        L_idx = findfirst(x -> qpts[1, 1, x][3] >= l_def, 1:size(qpts, 3))  # Find L index
        mode_slice = [modes[i, K_idx, L_idx] for i in 1:length(H_vals) if H_vals[i] in H_range]  # Slice along H
        q_values = H_range
    elseif typeof(k_def) <: Tuple
        # We're slicing along K (k_def is a range)
        K_vals = [qpts[1, j, 1][2] for j in 1:size(qpts, 2)]  # Extract K values from qpts
        K_range = filter(k -> k >= k_def[1] && k <= k_def[2], K_vals)  # Create K range within the slice range
        H_idx = findfirst(x -> qpts[x, 1, 1][1] >= h_def, 1:size(qpts, 1))  # Find H index
        L_idx = findfirst(x -> qpts[1, 1, x][3] >= l_def, 1:size(qpts, 3))  # Find L index
        mode_slice = [modes[H_idx, j, L_idx] for j in 1:length(K_vals) if K_vals[j] in K_range]  # Slice along K
        q_values = K_range
    elseif typeof(l_def) <: Tuple
        # We're slicing along L (l_def is a range)
        L_vals = [qpts[1, 1, k][3] for k in 1:size(qpts, 3)]  # Extract L values from qpts
        L_range = filter(l -> l >= l_def[1] && l <= l_def[2], L_vals)  # Create L range within the slice range
        H_idx = findfirst(x -> qpts[x, 1, 1][1] >= h_def, 1:size(qpts, 1))  # Find H index
        K_idx = findfirst(x -> qpts[1, x, 1][2] >= k_def, 1:size(qpts, 2))  # Find K index
        mode_slice = [modes[H_idx, K_idx, k] for k in 1:length(L_vals) if L_vals[k] in L_range]  # Slice along L
        q_values = L_range
    end

    return q_values, mode_slice
end

# Plot the slices for all available modes
function plot_sliced_modes(q_vals, mode_slices, title_str)
    p = plot(q_vals, [m[1] for m in mode_slices], label="Mode 1", xlabel="Momentum", ylabel="ωq", title=title_str)
    
    # Loop through remaining modes and plot them
    for mode_idx in 2:length(mode_slices[1])
        plot!(q_vals, [m[mode_idx] for m in mode_slices], label="Mode $mode_idx")
    end
    return p
end

# Helper function to generate a label based on the slice definition
function generate_label(slice_definition)
    h_def, k_def, l_def = slice_definition
    return "K=$(k_def), L=$(l_def)"
end

# Parameters for the momenta
H_vals = range(0, stop=2, length=101)
K_vals = range(0, stop=2, length=101)
L_vals = range(0, stop=2, length=101)

# Generate the qpts object containing all (H, K, L) tuples
qpts = generate_qpts(H_vals, K_vals, L_vals)

# Calculate the modes using qpts and package them into a single variable
modes = calculate_modes_from_qpts(qpts)

# Slice for K = 0, L = 0 along H
slice_definition_1 = [(0, 2.0), 0, 0]  # Slice along H
q_vals_1, modes_1_slice = oneD_mode_slice(modes, qpts, slice_definition_1)

# Slice for K = 0.5, L = 0.5 along H
slice_definition_2 = [(0, 2), 0.5, 0.5]  # Slice along H
q_vals_2, modes_2_slice = oneD_mode_slice(modes, qpts, slice_definition_2)

# Generate dynamic labels based on the slice definitions
label_1 = generate_label(slice_definition_1)
label_2 = generate_label(slice_definition_2)

# Create the plots
plot(
    plot_sliced_modes(q_vals_1, modes_1_slice, label_1),
    plot_sliced_modes(q_vals_2, modes_2_slice, label_2),
    layout=(2, 1)  # 2x1 layout for the plots
)
