function rand_sphere(args...)
    theta = acos.(1 .- 2 .* rand(args...))
    phi = 2π .* rand(args...)
    vcat(theta', phi')
end
