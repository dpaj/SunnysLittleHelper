using FFTW, Plots

# Function to create a 3D Gaussian kernel centered in q-space of the histogram
function create_gaussian_kernel(FWHM, qpts)
    # Convert FWHM to sigma
    sigma = FWHM / 2.35482

    # Initialize the kernel with the same size as qpts
    kernel = zeros(size(qpts))

    # Find the extrema (min, max) for H, K, L
    H_vals = [qpts[i, 1, 1][1] for i in 1:size(qpts, 1)]
    K_vals = [qpts[1, j, 1][2] for j in 1:size(qpts, 2)]
    L_vals = [qpts[1, 1, k][3] for k in 1:size(qpts, 3)]

    H_min, H_max = extrema(H_vals)
    K_min, K_max = extrema(K_vals)
    L_min, L_max = extrema(L_vals)

    # Calculate the midpoints for H, K, and L
    H_mid = (H_min + H_max) / 2
    K_mid = (K_min + K_max) / 2
    L_mid = (L_min + L_max) / 2

    # Loop over all points in qpts and calculate the Gaussian value based on q-vector distance from the midpoints
    for i in 1:size(qpts, 1)
        for j in 1:size(qpts, 2)
            for k in 1:size(qpts, 3)
                # Get the (H, K, L) values for this point in qpts
                H, K, L = qpts[i, j, k]
                
                # Calculate the squared distance in q-space relative to the midpoints
                r_squared = (H - H_mid)^2 + (K - K_mid)^2 + (L - L_mid)^2

                # Apply the Gaussian formula based on the distance from the numerical midpoints
                kernel[i, j, k] = exp(-r_squared / (2 * sigma^2))
            end
        end
    end
    
    # Normalize the kernel
    return kernel / sum(kernel)
end


# Function to perform real-to-complex FFT convolution of the 4D intensities in q-space
function convolve_4d_intensities_rfft(intensities_4d, FWHM, qpts)
    # Create the 3D Gaussian kernel for momentum broadening using qpts
    gaussian_kernel = create_gaussian_kernel(FWHM, qpts)

    # Initialize the convolved intensities array
    intensities_4d_convolved = similar(intensities_4d)

    # Perform the FFT convolution for each energy slice using rfft
    for e_idx in 1:size(intensities_4d, 4)
        # Get the energy slice from intensities_4d
        slice = intensities_4d[:, :, :, e_idx]

        # Perform real-to-complex FFT convolution
        fft_slice = rfft(slice)  # Forward rfft of the intensity slice
        fft_kernel = rfft(gaussian_kernel)  # Forward rfft of the Gaussian kernel

        # Multiply in Fourier space and then apply inverse real FFT (irfft)
        convolved_slice = irfft(fft_slice .* fft_kernel, size(slice, 1))  # Inverse FFT to real space

        # Store the convolved slice
        intensities_4d_convolved[:, :, :, e_idx] = convolved_slice
    end

    return intensities_4d_convolved
end

# Helper function to find the closest index to a given value
function find_closest_idx(val, arr)
    return argmin(abs.(arr .- val))  # Find the index of the closest value
end

# Parameters for the momenta
H_vals = range(0, stop=2, length=101)
K_vals = range(0, stop=2, length=101)
L_vals = range(0, stop=2, length=101)

# Generate the qpts object containing all (H, K, L) tuples
qpts = Array{Tuple{Float64, Float64, Float64}}(undef, length(H_vals), length(K_vals), length(L_vals))
for i in 1:length(H_vals)
    for j in 1:length(K_vals)
        for k in 1:length(L_vals)
            qpts[i, j, k] = (H_vals[i], K_vals[j], L_vals[k])
        end
    end
end


# Set the FWHM for momentum broadening (same for H, K, and L)
momentum_FWHM = 0.1

@time begin
    # Perform the convolution on the 4D intensities (HKL + E)
    intensities_4d_convolved = convolve_4d_intensities_rfft(intensities_4d, momentum_FWHM, qpts)
end

# Define energy points for plotting
energies = range(0, stop=30, length=num_energy_points)

# Example: User inputs qpt values
qpt_input = [0.5, 0.5, 0.25]  # Example user input for H, K, L

# Find the closest indices for the given qpt values
H_idx = find_closest_idx(qpt_input[1], H_vals)
K_idx = find_closest_idx(qpt_input[2], K_vals)
L_idx = find_closest_idx(qpt_input[3], L_vals)

# Extract the intensity at the found indices from the convolved intensities
intensity_at_point_convolved = intensities_4d_convolved[H_idx, K_idx, L_idx, :]

# Plot the convolved intensities for that momentum point as a function of energy
plot(energies, intensity_at_point_convolved, label="Intensity (Convolved)", xlabel="Energy", ylabel="Intensity",
     title="Intensity at (H=$(H_vals[H_idx]), K=$(K_vals[K_idx]), L=$(L_vals[L_idx])) with Q-space broadening")
