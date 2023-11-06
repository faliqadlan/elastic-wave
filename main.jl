using Plots

# Parameters
c = 850     # Speed of wave propagation (in m/s)
f = 1.0        # Central frequency of the source (in Hz)
λ = c / f      # Wavelength

# Adjust the spatial and temporal parameters to satisfy the CFL condition
dx = 0.01 * λ  # Spatial step size (about 1/10 of a wavelength)
dt = dx / c    # Time step size based on the CFL condition

t_max = 10.0    # Maximum simulation time (in seconds)
x_max = 10.0    # Maximum spatial coordinate (in meters)

nt = round(Int, t_max / dt) # Calculate the number of time steps based on t_max
nx = round(Int, x_max / dx) # Calculate the number of spatial steps based on x_max

ρ = 2.5        # Density of the medium (in kg/m^3)
μ = ρ * c^2    # Shear modulus of the medium (in Pa)

u = zeros(nx)       # Displacement at current time step
u_prev = zeros(nx)  # Displacement at the previous time step
u_next = zeros(nx)  # Displacement at the next time step

u_time = zeros(nt, nx)  # Displacement at all time steps

# Source term (simple Ricker wavelet)
function ricker_wavelet(t, f; tau=2.0)
    A = (1 - 2 * π^2 * f^2 * (t - tau)^2)
    wavelet = A * exp(-π^2 * f^2 * (t - tau)^2)
    return wavelet
end

source_loc = nx ÷ 2

# Time-stepping loop
for it = 1:nt
    for ix = 2:nx-1
        u_next[ix] = 2u[ix] - u_prev[ix] + (c^2 * dt^2 / dx^2) * (u[ix+1] - 2u[ix] + u[ix-1])
    end
    u_next[source_loc] += ricker_wavelet(it * dt, f)

    u_prev, u, u_next = u, u_next, u_prev

    u_time[it, :] = u
end

# max_u = maximum(u_time)
# min_u = minimum(u_time)

# for it = 1:nt
#     IJulia.clear_output(true)
#     p = plot(u_time[it, :], legend=false, title="time: $(it*dt) s", ylims=(min_u, max_u))
#     display(p)
#     sleep(0.1)
# end

# anim = @animate for it = 1:nt
#     p = plot(u_time[it, :], legend=false, title="time: $(it*dt) s", ylims=(min_u, max_u))
# end

# gif(anim, "elastic_wave.gif", fps=20)