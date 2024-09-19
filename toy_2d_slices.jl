using Plots

# Function to extract a 2D slice based on user input for momentum and energy
function extract_2d_slice(intensities_4d, H_vals, K_vals, L_vals, energies, slice_definition)
    h_def, k_def, l_def, e_def = slice_definition  # The user-provided slice definition

    # Determine which dimensions are 2D ranges and which are fixed
    if typeof(h_def) <: Tuple && typeof(e_def) <: Tuple
        # Slice along H and E, fix K and L
        K_idx = find_closest_idx(k_def, K_vals)
        L_idx = find_closest_idx(l_def, L_vals)

        # Map H and E to their nearest indices
        H_start_idx = find_closest_idx(h_def[1], H_vals)
        H_end_idx = find_closest_idx(h_def[2], H_vals)
        E_start_idx = find_closest_idx(e_def[1], energies)
        E_end_idx = find_closest_idx(e_def[2], energies)

        # Extract the 2D slice (H vs E) at fixed K and L
        slice = intensities_4d[H_start_idx:H_end_idx, K_idx, L_idx, E_start_idx:E_end_idx]
        x_vals = H_vals[H_start_idx:H_end_idx]  # H values for the x-axis
        y_vals = energies[E_start_idx:E_end_idx]  # Energy values for the y-axis

        return x_vals, y_vals, slice

    elseif typeof(h_def) <: Tuple && typeof(k_def) <: Tuple
        # Slice along H and K, fix L and E
        L_idx = find_closest_idx(l_def, L_vals)
        E_idx = find_closest_idx(e_def, energies)

        # Map H and K to their nearest indices
        H_start_idx = find_closest_idx(h_def[1], H_vals)
        H_end_idx = find_closest_idx(h_def[2], H_vals)
        K_start_idx = find_closest_idx(k_def[1], K_vals)
        K_end_idx = find_closest_idx(k_def[2], K_vals)

        # Extract the 2D slice (H vs K) at fixed L and E
        slice = intensities_4d[H_start_idx:H_end_idx, K_start_idx:K_end_idx, L_idx, E_idx]
        x_vals = H_vals[H_start_idx:H_end_idx]  # H values for the x-axis
        y_vals = K_vals[K_start_idx:K_end_idx]  # K values for the y-axis

        return x_vals, y_vals, slice

    else
        error("Invalid slice definition. Provide two ranges and two fixed values.")
    end
end

# First user-provided slice definition (for example, slice along H and E with fixed K and L)
slice_definition1 = [(0.5, 1.5), 0.75, 1.0, (4.0, 25.0)]  # Slice along H and E with fixed K and L

# Second user-provided slice definition (for example, slice along H and K with fixed L and E)
slice_definition2 = [(0.0, 2.0), (0.0, 2.0), 1.0, 10.0]  # Slice along H and K with fixed L and E

# Extract the 2D slice from original data for first slice definition
x_vals1, y_vals1, slice1 = extract_2d_slice(intensities_4d, H_vals, K_vals, L_vals, energies, slice_definition1)

# Extract the 2D slice from convolved data for first slice definition
x_vals2, y_vals2, slice2 = extract_2d_slice(intensities_4d_convolved, H_vals, K_vals, L_vals, energies, slice_definition1)

# Extract the 2D slice from original data for second slice definition
x_vals3, y_vals3, slice3 = extract_2d_slice(intensities_4d, H_vals, K_vals, L_vals, energies, slice_definition2)

# Extract the 2D slice from convolved data for second slice definition
x_vals4, y_vals4, slice4 = extract_2d_slice(intensities_4d_convolved, H_vals, K_vals, L_vals, energies, slice_definition2)

# Create the compound plot with 4 subplots (2x2 layout)
p1 = heatmap(x_vals1, y_vals1, slice1', xlabel="H", ylabel="Energy", title="Original: Slice 1 (H vs E)")
p2 = heatmap(x_vals2, y_vals2, slice2', xlabel="H", ylabel="Energy", title="Convolved: Slice 1 (H vs E)")
p3 = heatmap(x_vals3, y_vals3, slice3', xlabel="H", ylabel="K", title="Original: Slice 2 (H vs K)")
p4 = heatmap(x_vals4, y_vals4, slice4', xlabel="H", ylabel="K", title="Convolved: Slice 2 (H vs K)")

# Display the compound plot
plot(p1, p2, p3, p4, layout=(2, 2), size=(1000, 800))
