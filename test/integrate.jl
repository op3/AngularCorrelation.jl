using FastGaussQuadrature

"""
Integrate f(a, b) over a, b ∈ [-1, 1]
using Gauss quadrature with n-point nodes.
"""
function Quad2D(f, n)
    x, w = gausslegendre(n)
    y = x
    sum = 0
    for i = 1:n, j = 1:n
        sum = sum + f(x[i], y[j]) * w[i] * w[j]
    end
    return sum
end

"""
Integrate f(theta, phi) over theta ∈ [0, π], phi ∈ [0, 2π]
using Gauss quadrature with n-point nodes.
"""
function QuadSphere(f, n)
    Quad2D(
        function (a, b)
            theta = (a + 1) * π / 2  # [0, π]
            phi = (b + 1) * π  # [0, 2π]
            return sin(theta) * π * π / 2 * f(theta, phi)
        end, n)
end

"""
Integrate f(a, b, c, d) over a, b, c, d ∈ [-1, 1]
using Gauss quadrature with n-point nodes.
"""
function Quad4D(f, n)
    x0, w = gausslegendre(n)
    x1 = x0
    x2 = x0
    x3 = x0
    sum = 0
    for i = 1:n, j = 1:n, k = 1:n, l = 1:n
        sum = sum + f(x0[i], x1[j], x2[k], x3[l]) * w[i] * w[j] * w[k] * w[l]
    end
    return sum
end

"""
Integrate f(theta1, phi1, theta2, phi2) over
theta1, theta2 ∈ [0, π], phi1, phi2 ∈ [0, 2π]
using Gauss quadrature with n-point nodes.
"""
function QuadSphere2(f, n)
    Quad4D(
        function (a, b, c, d)
            theta1 = (a + 1) * π / 2  # [0, π]
            phi1 = (b + 1) * π  # [0, 2π]
            theta2 = (c + 1) * π / 2  # [0, π]
            phi2 = (d + 1) * π  # [0, 2π]
            return sin(theta1) * π * π / 2 * sin(theta2) * π * π / 2 * f(theta1, phi1, theta2, phi2)
        end, n)
end

