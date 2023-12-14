using LinearAlgebra
using Plots
using SparseArrays

# Basic parameters
vs = 100.0  # Wave velocity [m/s]
f = 0.5  # Central frequency of the source (in Hz)
λ = vs / f  # Wavelength

# Spatial and temporal parameters based on CFL condition
ep = 0.3  # Stability limit
dx = 0.1 * λ  # Spatial step size
dx_min = dx
dt = ep * dx / vs  # Time step size based on CFL

t_max = 0.5  # Maximum simulation time (in seconds)
x_max = 100.0  # Maximum spatial coordinate (in meters)

# Print information
println("Spatial step size: $dx")
println("Time step size: $dt")

# Calculate number of steps
nt = round(Int, t_max / dt)
nx = round(Int, x_max / dx)
ny = nx  # Assuming square domain for 2D
println("Number of time steps: $nt")
println("Number of spatial steps: $nx, $ny")

# Initialize spatial coordinates
elem = zeros(nx, ny)

for i in 1:nx
    for j in 1:ny
        elem[i, j] = (i - 1) * nx + j
    end
end


# parameters physical
ro0 = 1500    # Density [kg/m^3]
ro = elem * 0 .+ ro0          # initialize density array
young_mod = elem * 0 .+ ro .* vs .^ 2   # calculate young modulus from density and velocity

# Time axis
t = range(dt, stop=nt * dt, length=nt)

# Displacement
u = zeros(nx, ny)
uold = zeros(nx, ny)
unew = zeros(nx, ny)

# Source term (simple Ricker wavelet)
function ricker_wavelet(t, f, tau=2.0)
    # A = (1 - 2 * π^2 * f^2 * (t - tau)^2)
    # wavelet = A * exp(-π^2 * f^2 * (t - tau)^2)
    wavelet = sin(2 * π * f * (t - tau)) * exp(-π^2 * f^2 * (t - tau)^2)
    return wavelet
end

function bilinear_basis(xi, yi)
    return [(1 - xi) * (1 - yi) / 4, (1 + xi) * (1 - yi) / 4, (1 + xi) * (1 + yi) / 4, (1 - xi) * (1 + yi) / 4]
end

function jacobian(x1, x2, y1, y2)
    # Calculate basis function values at each point
    b11 = (1 - x1 / (x2 - x1)) * (1 - y1 / (y2 - y1)) / 4
    b12 = (1 + x1 / (x2 - x1)) * (1 - y1 / (y2 - y1)) / 4
    b21 = (1 + x1 / (x2 - x1)) * (1 + y1 / (y2 - y1)) / 4
    b22 = (1 - x1 / (x2 - x1)) * (1 + y1 / (y2 - y1)) / 4

    # Define Jacobian matrix
    J = zeros(2, 2)
    J[1, 1] = (b22 - b12) / (x2 - x1)
    J[1, 2] = (b21 - b11) / (y2 - y1)
    J[2, 1] = (b12 - b22) / (x2 - x1)
    J[2, 2] = (b11 - b21) / (y2 - y1)

    return J
end

function gauss_legendre_quadrature()
    # Define quadrature points and weights for 2d
    points = [
        (-sqrt((3 / 7) - ((2 / 7) * sqrt(6 / 5))), -sqrt((3 / 7) - ((2 / 7) * sqrt(6 / 5)))),
        (sqrt((3 / 7) - ((2 / 7) * sqrt(6 / 5))), -sqrt((3 / 7) - ((2 / 7) * sqrt(6 / 5)))),
        (-sqrt((3 / 7) + ((2 / 7) * sqrt(6 / 5))), sqrt((3 / 7) + ((2 / 7) * sqrt(6 / 5)))),
        (sqrt((3 / 7) + ((2 / 7) * sqrt(6 / 5))), sqrt((3 / 7) + ((2 / 7) * sqrt(6 / 5))))
    ]

    weights = [(18 + sqrt(30)) / 36, (18 + sqrt(30)) / 36, (18 - sqrt(30)) / 36, (18 - sqrt(30)) / 36]

    return points, weights
end


function mass_matrix(x1, x2, y1, y2, ro, A=1.0)
    points, weights = gauss_legendre_quadrature()

    M = spzeros(4, 4)  # Assuming bilinear elements, hence 4x4 mass matrix

    for k in 1:4
        xi, yi = points[k]
        w = weights[k]

        basis_values = bilinear_basis(xi, yi)

        for i in 1:4
            for j in 1:4
                M[i, j] += w * basis_values[i] * basis_values[j] * ro * A * det(jacobian(x1, x2, y1, y2))
            end
        end
    end

    return M
end

function stiffness_matrix(x1, x2, y1, y2, young_mod)
    points, weights = gauss_legendre_quadrature()

    num_nodes = 4
    K = spzeros(num_nodes, num_nodes)  # Assuming bilinear elements, hence 4x4 stiffness matrix

    basis_derivatives = [-1 / 2, 1 / 2]  # Derivatives of linear basis functions

    for k in 1:4
        xi, yi = points[k]
        w = weights[k]

        for i in 1:num_nodes
            for j in 1:num_nodes
                dN_dx = basis_derivatives[i] / jacobian(x1, x2, y1, y2)

                K[i, j] += w * young_mod * (dN_dx * dN_dx') * jacobian(x1, x2, y1, y2)
            end
        end
    end

    return K
end

function global_matrix(elem, ro, young_mod)
    num_nodes = size(elem, 1)
    num_elements = (size(elem, 1) - 1) * (size(elem, 2) - 1)  # Assuming a regular grid

    M_global = spzeros(num_nodes, num_nodes)
    K_global = spzeros(num_nodes, num_nodes)

    for i in 1:size(elem, 1)-1
        for j in 1:size(elem, 2)-1
            x1, y1 = elem[i, j], elem[i, j+1]
            x2, y2 = elem[i+1, j], elem[i+1, j+1]

            M_e = mass_matrix(x1, x2, y1, y2, ro[i, j])
            K_e = stiffness_matrix(x1, x2, y1, y2, young_mod[i, j])

            indices = [x1, y1, x2, y2]  # Nodes associated with the element

            M_global[indices, indices] += M_e
            K_global[indices, indices] += K_e
        end
    end

    return M_global, K_global
end


# Invert M using pseudoinverse
mass_matrices, stiffness_matrices = global_matrix(elem, ro, young_mod)
mass_matrices = transpose(mass_matrices)
stiffness_matrices = transpose(stiffness_matrices)

# Plot mass matrix and stiffness matrix in one plot
plot(
    heatmap(mass_matrices, title="Mass Matrix Inv", color=:viridis, legend=false, yflip=true),
    heatmap(stiffness_matrices, title="Stiffness Matrix K", color=:viridis, legend=false, yflip=true),
    layout=(2, 1), size=(400, 800)
)
