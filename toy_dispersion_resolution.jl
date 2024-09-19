using Plots

# Dispersion relation with Delta
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
    # We will now store 3 modes for demonstration (can be any number of modes)
    modes = Array{NTuple{3, Float64}}(undef, size(qpts)...)  # Each element will store multiple modes
    
    for i in 1:size(qpts, 1)
        for j in 1:size(qpts, 2)
            for k in 1:size(qpts, 3)
                H, K, L = qpts[i, j, k]
                mode1 = dispersion_relation(5, 1, H, K, L, 1e-3)  # Mode 1: S = 1, J = 5, Δ ≈ 0
                mode2 = dispersion_relation(2, 1, H, K, L, 5)  # Mode 2: S = 1, J = 2, Δ = 5
                mode3 = dispersion_relation(3, 1, H, K, L, 2)  # Mode 3: S = 1, J = 3, Δ = 2
                modes[i, j, k] = (mode1, mode2, mode3)  # Store all modes as a tuple
            end
        end
    end
    return modes
end

# Gaussian resolution function depending on energy
function resolution(energy)
    #my_fwhm =  3.8775e-05 * energy^3 + 0.0018964 * energy^2 - 0.078275 * energy + 0.72305 #12 meV Ei, CNCS
    my_fwhm =  +7.931e-06 * energy^3 +0.0011244 * energy^2 -0.12422 * energy +2.7828 #30 meV Ei, CNCS
    #my_fwhm =  1.0 #fake constant width
    return my_fwhm/2.35482
end

# Create a 4D intensities object (H, K, L, energy) and fill it based on the modes and resolution
function calculate_4d_intensities_first_try(modes, qpts, energy_range, num_energy_points)
    intensities = zeros(size(qpts)..., num_energy_points)  # 4D array for intensities
    energies = range(energy_range[1], energy_range[2], length=num_energy_points)  # Energy points
    
    for i in 1:size(modes, 1)
        for j in 1:size(modes, 2)
            for k in 1:size(modes, 3)
                for mode in modes[i, j, k]  # Iterate over each mode at a given (H, K, L)
                    # Get the width from the resolution function
                    width = resolution(mode)
                    # Compute Gaussian contributions for this mode, the 1/energy dependence of intensity is typical for an AFM magnon
                    for e_idx in 1:num_energy_points
                        energy = energies[e_idx]
                        intensity = (1 / mode) * exp(-((energy - mode)^2) / (2 * width^2))  # Mode intensity with Gaussian, AFM
                        #intensity =  exp(-((energy - mode)^2) / (2 * width^2))  # Mode intensity with Gaussian, FM
                        intensities[i, j, k, e_idx] += intensity  # Accumulate the intensity for each energy point
                    end
                end
            end
        end
    end
    
    return energies, intensities
end

using Base.Threads  # To enable multi-threading

# Optimized function to calculate the 4D intensities
function calculate_4d_intensities_trying_speedup(modes, qpts, energy_range, num_energy_points)
    intensities = zeros(size(qpts)..., num_energy_points)  # 4D array for intensities
    energies = range(energy_range[1], energy_range[2], length=num_energy_points)  # Energy points
    
    # Parallelize over the outer momentum loops
    @threads for i in 1:size(modes, 1)
        for j in 1:size(modes, 2)
            for k in 1:size(modes, 3)
                for mode in modes[i, j, k]  # Iterate over each mode at a given (H, K, L)
                    # Pre-calculate the resolution width for this mode
                    width = resolution(mode)
                    if width == 0  # Avoid division by zero
                        continue
                    end
                    # Pre-compute the intensity factor for this mode, the 1/energy dependence is typical for an AFM magnon
                    mode_intensity_factor = 1 / mode
                    #mode_intensity_factor = 1
                    # Pre-compute the Gaussian contribution for all energy points (vectorized)
                    energy_diffs = energies .- mode
                    gaussian_contribution = exp.(-((energy_diffs .^ 2) / (2 * width^2)))
                    intensity_contribution = mode_intensity_factor * gaussian_contribution
                    
                    # Add the intensity contribution for all energy points at once
                    intensities[i, j, k, :] .+= intensity_contribution
                end
            end
        end
    end
    
    return energies, intensities
end


# Perform 1D slicing along H, K, or L
function oneD_mode_slice(modes, qpts, slice_definition)
    h_def, k_def, l_def = slice_definition
    
    # Determine which dimension is a range and which are fixed
    if typeof(h_def) <: Tuple
        H_vals = [qpts[i, 1, 1][1] for i in 1:size(qpts, 1)]  # Extract H values from qpts
        H_range = filter(h -> h >= h_def[1] && h <= h_def[2], H_vals)  # Create H range within the slice range
        K_idx = findfirst(x -> qpts[1, x, 1][2] >= k_def, 1:size(qpts, 2))  # Find K index
        L_idx = findfirst(x -> qpts[1, 1, x][3] >= l_def, 1:size(qpts, 3))  # Find L index
        mode_slice = [modes[i, K_idx, L_idx] for i in 1:length(H_vals) if H_vals[i] in H_range]  # Slice along H
        q_values = H_range
    elseif typeof(k_def) <: Tuple
        K_vals = [qpts[1, j, 1][2] for j in 1:size(qpts, 2)]  # Extract K values from qpts
        K_range = filter(k -> k >= k_def[1] && k <= k_def[2], K_vals)  # Create K range within the slice range
        H_idx = findfirst(x -> qpts[x, 1, 1][1] >= h_def, 1:size(qpts, 1))  # Find H index
        L_idx = findfirst(x -> qpts[1, 1, x][3] >= l_def, 1:size(qpts, 3))  # Find L index
        mode_slice = [modes[H_idx, j, L_idx] for j in 1:length(K_vals) if K_vals[j] in K_range]  # Slice along K
        q_values = K_range
    elseif typeof(l_def) <: Tuple
        L_vals = [qpts[1, 1, k][3] for k in 1:size(qpts, 3)]  # Extract L values from qpts
        L_range = filter(l -> l >= l_def[1] && l <= l_def[2], L_vals)  # Create L range within the slice range
        H_idx = findfirst(x -> qpts[x, 1, 1][1] >= h_def, 1:size(qpts, 1))  # Find H index
        K_idx = findfirst(x -> qpts[1, x, 1][2] >= k_def, 1:size(qpts, 2))  # Find K index
        mode_slice = [modes[H_idx, K_idx, k] for k in 1:length(L_vals) if L_vals[k] in L_range]  # Slice along L
        q_values = L_range
    end

    return q_values, mode_slice
end



# Plot the slices dynamically for all available modes
function plot_sliced_modes(q_vals, mode_slices, title_str)
    p = plot(q_vals, [m[1] for m in mode_slices], label="Mode 1", xlabel="Momentum", ylabel="ωq", title=title_str)
    
    # Loop through remaining modes and plot them
    for mode_idx in 2:length(mode_slices[1])
        plot!(q_vals, [m[mode_idx] for m in mode_slices], label="Mode $mode_idx")
    end
    return p
end

# Function to find the index of the closest value in an array
function find_closest_idx(val, arr)
    return argmin(abs.(arr .- val))  # Find the index of the closest value
end

# Parameters for the momenta
H_vals = range(0, stop=2, length=101)
K_vals = range(0, stop=2, length=101)
L_vals = range(0, stop=2, length=101)

# Generate the qpts object containing all (H, K, L) tuples
qpts = generate_qpts(H_vals, K_vals, L_vals)

# Calculate the modes using qpts and package them into a single variable
modes = calculate_modes_from_qpts(qpts)

# Energy range and number of points for intensities
energy_range = (0.0, 30.0)
num_energy_points = 301

# Fill the 4-d object HKLE with intensities broadened by the instrumental resolution, this is expensive (for now?)...
@time begin
#energies, intensities_4d = calculate_4d_intensities_first_try(modes, qpts, energy_range, num_energy_points)
energies, intensities_4d = calculate_4d_intensities_trying_speedup(modes, qpts, energy_range, num_energy_points)
end

# User inputs qpt values
qpt_input = [0.5, 0.5, 0.25]  # Example user input for H, K, L, play around with this

# Find the closest indices for the given qpt values
H_idx = find_closest_idx(qpt_input[1], H_vals)
K_idx = find_closest_idx(qpt_input[2], K_vals)
L_idx = find_closest_idx(qpt_input[3], L_vals)

# Extract the intensity at the found indices
intensity_at_point = intensities_4d[H_idx, K_idx, L_idx, :]

# Plot the intensities for that momentum point as a function of energy
plot(energies, intensity_at_point, label="Intensity", xlabel="Energy", ylabel="Intensity",
     title="Intensity at (H=$(H_vals[H_idx]), K=$(K_vals[K_idx]), L=$(L_vals[L_idx]))")

