using Plots

# Parameters
nx = 20000     # Number of spatial points in the domain
dx = 1.0       # Spatial step size (in meters)
dt = 0.05    # Time step size (in seconds)
nt = 100      # Number of time steps to simulate
ρ = 2.5        # Density of the medium (in kg/m^3)
μ = 2e6        # Shear modulus of the medium (in Pa)
c = sqrt(μ / ρ)  # Speed of wave propagation (in m/s)
f = 1.0        # Central frequency of the source (in Hz)

# Source term (simple Ricker wavelet)
function ricker_wavelet(t, f; tau=2.0)
    A = (1 - 2 * π^2 * f^2 * (t - tau)^2)
    wavelet = A * exp(-π^2 * f^2 * (t - tau)^2)
    return wavelet
end

source_loc = nx ÷ 2

# Create an empty plot for later use
p = plot(legend=false, title="time: 0.0 s", xlims=(0, nx))

global u = zeros(nx)       # Displacement at current time step
global u_prev = zeros(nx)  # Displacement at previous time step
global u_next = zeros(nx)  # Displacement at next time step

# Time-stepping loop
anim = @animate for it = 1:nt

    for ix = 2:nx-1
        u_next[ix] = 2u[ix] - u_prev[ix] + (c^2 * dt^2 / dx^2) * (u[ix+1] - 2u[ix] + u[ix-1])
    end
    u_next[source_loc] += ricker_wavelet(it * dt, f)

    u_prev, u, u_next = u, u_next, u_prev

    # Update the plot at every 10 time steps
    if it % 10 == 0
        plot!(p, u, legend=false, title="time: $(it*dt)")
    end
end

gif(anim, "wave_propagation.gif", fps=15)