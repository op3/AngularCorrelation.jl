using AngularCorrelation
using FastGaussQuadrature
using Test

include("integrate.jl")

@testset "AngularCorrelation.jl" begin

    @test E1() ≡ Transition(1, EMCharacter.electric, 2, EMCharacter.magnetic, 0.0)
    @test M1() ≡ Transition(1, EMCharacter.magnetic, 2, EMCharacter.electric, 0.0)
    @test Dipole() ≡ Transition(1, EMCharacter.unknown, 2, EMCharacter.unknown, 0.0)
    @test E2() ≡ Transition(2, EMCharacter.electric, 3, EMCharacter.magnetic, 0.0)
    @test M2() ≡ Transition(2, EMCharacter.magnetic, 3, EMCharacter.electric, 0.0)
    @test Quadrupole() ≡ Transition(2, EMCharacter.unknown, 3, EMCharacter.unknown, 0.0)
    @test E1(3.3) ≡ Transition(1, EMCharacter.electric, 2, EMCharacter.magnetic, 3.3)
    @test M1(-3.3) ≡ Transition(1, EMCharacter.magnetic, 2, EMCharacter.electric, -3.3)
    @test Dipole(2.2) ≡ Transition(1, EMCharacter.unknown, 2, EMCharacter.unknown, 2.2)
    @test E2(-1.0) ≡ Transition(2, EMCharacter.electric, 3, EMCharacter.magnetic, -1.0)
    @test M2(1.0) ≡ Transition(2, EMCharacter.magnetic, 3, EMCharacter.electric, 1.0)
    @test Quadrupole(2.2) ≡ Transition(2, EMCharacter.unknown, 3, EMCharacter.unknown, 2.2)

    @test Int(EMCharacter.electric) ≡ (-1)^0
    @test Int(EMCharacter.magnetic) ≡ (-1)^1

    theta1 = π / 3
    phi1 = π / 3
    theta2 = π / 3
    phi2 = 2 * π / 3

    # 0⁺ → 1⁺ → 2 → 0
    @test Wcorr(
        theta1, phi1, theta2, phi2,
        State(0), M1(),
        State(1), Dipole(),
        State(2), Quadrupole(),
        State(0)) ≈ 0.0744873 rtol = 1e-5

    # 0⁺ → 1⁻ → 2 → 0
    @test Wcorr(
        theta1, phi1, theta2, phi2,
        State(0), E1(),
        State(1), Dipole(),
        State(2), Quadrupole(),
        State(0)) ≈ 0.0975679 rtol = 1e-5

    # 0⁺ → 2⁺ → 2 → 0
    @test Wcorr(
        theta1, phi1, theta2, phi2,
        State(0), E2(),
        State(2), Dipole(),
        State(2), Quadrupole(),
        State(0)) ≈ 0.0367191 rtol = 1e-5

    # 0⁺ → 1⁺ → 2 → 2
    @test Wcorr(
        theta1, phi1, theta2, phi2,
        State(0), M1(),
        State(1), Dipole(),
        State(2), Dipole(),
        State(2)) ≈ 0.0889651 rtol = 1e-5

    # 0⁺ → 1⁻ → 2 → 2
    @test Wcorr(
        theta1, phi1, theta2, phi2,
        State(0), E1(),
        State(1), Dipole(),
        State(2), Dipole(),
        State(2)) ≈ 0.077215 rtol = 1e-5

    # 0⁺ → 2⁺ → 2 → 2
    @test Wcorr(
        theta1, phi1, theta2, phi2,
        State(0), E2(),
        State(2), Dipole(),
        State(2), Dipole(),
        State(2)) ≈ 0.0562638 rtol = 1e-5

    @test QuadSphere2(
        function (theta1, phi1, theta2, phi2)
            return Wcorr(
                theta1, phi1, theta2, phi2,
                State(0), E2(),
                State(2), Dipole(),
                State(2), Dipole(),
                State(2))
        end, 9) ≈ 4 * π rtol = 1e-6

    @test QuadSphere2(
        function (theta1, phi1, theta2, phi2)
            return Wcorr(
                theta1, phi1, theta2, phi2,
                State(0), E2(),
                State(2), Dipole(),
                State(2), Dipole(),
                State(2), Dipole(),
                State(2))
        end, 9) ≈ 4 * π rtol = 1e-6

    @test QuadSphere(
        function (theta, phi)
            return W(
                theta, phi,
                State(0), E2(),
                State(2), Quadrupole(),
                State(0))
        end, 9) ≈ 4 * π rtol = 1e-5

    @test QuadSphere(
        function (theta, phi)
            return W(
                theta, phi,
                State(0), E2(),
                State(2), Dipole(),
                State(2), Quadrupole(),
                State(0))
        end, 9) ≈ 4 * π rtol = 1e-5

    @test W(π / 3, π / 3, State(0), M1(), State(1), Dipole(), State(0)) ≈ 0.65625
    @test W(π / 3, π / 3, State(0), E1(), State(1), Dipole(), State(0)) ≈ 1.21875
    @test W(π / 3, π / 3, State(0), M1(), State(1), Dipole(), State(2)) ≈ 0.965625
    @test W(π / 3, π / 3, State(0), E1(), State(1), Dipole(), State(2)) ≈ 1.021875
    @test W(π / 3, π / 3, State(0), M1(), State(1), Dipole(), State(2), Dipole(), State(2)) ≈ 1.1203125
    @test W(π / 3, π / 3, State(0), E1(), State(1), Dipole(), State(2), Dipole(), State(2)) ≈ 0.9234375
    @test W(π / 3, π / 3, State(0), E1(), State(1), Dipole(), State(2), Dipole(), State(2), Dipole(), State(2)) ≈ 0.96171875
    @test W(π / 3, π / 3, State(0), E1(), State(1), Dipole(), State(2), Dipole(), State(2), Quadrupole(), State(0)) ≈ 0.9453125
    @test W(π / 3, π / 3, State(0), E1(), State(1), Dipole(), State(2), Quadrupole(), State(0), Quadrupole(), State(2)) ≈ 1.0
    @test W(π / 3, π / 3, State(0), E1(), State(1), Dipole(1.0), State(2)) ≈ 1.21237 rtol = 1e-5
    @test W(π / 3, π / 3, State(0), E1(), State(1), Dipole(100.0), State(2)) ≈ 1.1123 rtol = 1e-6
    @test W(π / 3, π / 3, State(0), M1(), State(1), Dipole(100.0), State(2)) ≈ 0.823527 rtol = 1e-6
    @test W(π / 3, π / 3, State(0), M1(), State(1), Dipole(), State(2), Quadrupole(), State(0)) ≈ 1.171875
    @test W(π / 3, π / 3, State(0), M1(), State(1), Dipole(100.0), State(2), Quadrupole(), State(0)) ≈ 0.82815937
    @test W(π / 3, π / 3, State(0), M1(), State(1), Dipole(100.0), State(2), Dipole(), State(3)) ≈ 1.034368126
    @test W(π / 3, π / 3, State(0), M1(), State(1), Dipole(100.0), State(2), Dipole(100.0), State(3)) ≈ 1.101952726

    @test QuadSphere(
        function (theta, phi)
            return Wcorr(
                theta, phi,
                π / 5, π / 3,
                State(0), E1(),
                State(1), Dipole(),
                State(2), Quadrupole(),
                State(0))
        end, 9) ≈ W(
        π / 5, π / 3,
        State(0), E1(),
        State(1), Dipole(),
        State(2), Quadrupole(),
        State(0))

    @test QuadSphere(
        function (theta, phi)
            return Wcorr(
                theta, phi,
                π / 5, π / 3,
                State(0), E1(),
                State(1), Dipole(),
                State(2), Dipole(),
                State(2), Quadrupole(),
                State(0))
        end, 9) ≈ W(
        π / 5, π / 3,
        State(0), E1(),
        State(1), Dipole(),
        State(2), Dipole(),
        State(2), Quadrupole(),
        State(0))

    @test QuadSphere(
        function (theta, phi)
            return Wcorr(
                theta, phi,
                π / 5, π / 3, State(0), M1(),
                State(1), Dipole(),
                State(2), Quadrupole(),
                State(0))
        end, 9) ≈ W(
        π / 5, π / 3,
        State(0), M1(),
        State(1), Dipole(),
        State(2), Quadrupole(),
        State(0)) rtol = 1e-6

    @test QuadSphere(
        function (theta, phi)
            return Wcorr(
                theta, phi, π / 5, π / 5,
                State(0), M1(),
                State(1), Dipole(),
                State(2), Dipole(),
                State(2), Quadrupole(),
                State(0))
        end, 9) ≈ W(
        π / 5, π / 5,
        State(0), M1(),
        State(1), Dipole(),
        State(2), Dipole(),
        State(2), Quadrupole(),
        State(0))

    @test QuadSphere(
        function (theta, phi)
            return Wcorr(
                theta, phi, π / 3, π / 3,
                State(0), E1(), State(1), Dipole(100.0), State(2), Dipole(), State(3))
        end, 11) ≈ W(π / 3, π / 3, State(0), E1(), State(1), Dipole(100.0), State(2), Dipole(), State(3)) rtol = 1e-6

    @test QuadSphere(
        function (theta, phi)
            return Wcorr(
                theta, phi, π / 3, π / 3,
                State(0), E1(), State(1), Dipole(100.0), State(2), Dipole(100.0), State(3))
        end, 11) ≈ W(π / 3, π / 3, State(0), E1(), State(1), Dipole(100.0), State(2), Dipole(100.0), State(3)) rtol = 1e-5

    coeff = W_coeff(State(0), E1(), State(1), Dipole(), State(0))
    @test size(W_sample(10, coeff)) ≡ (2, 10)

    coeff = Wcorr_coeff(State(0), E1(), State(1), Dipole(), State(2), Quadrupole(), State(0))
    @test size(Wcorr_sample(100, coeff)) ≡ (4, 100)
    @test size(Wcorr_sample_mt(100, 2, coeff)) ≡ (4, 100)
end
