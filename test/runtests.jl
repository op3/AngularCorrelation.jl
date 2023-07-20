using AngularCorrelation
using Test

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
    State(0)) ≈ 0.0744873 rtol = 0.00001

# 0⁺ → 1⁻ → 2 → 0
@test Wcorr(
    theta1, phi1, theta2, phi2,
    State(0), E1(),
    State(1), Dipole(),
    State(2), Quadrupole(),
    State(0)) ≈ 0.0975679 rtol = 0.00001

# 0⁺ → 2⁺ → 2 → 0
@test Wcorr(
    theta1, phi1, theta2, phi2,
    State(0), E2(),
    State(2), Dipole(),
    State(2), Quadrupole(),
    State(0)) ≈ 0.0367191 rtol = 0.00001

# 0⁺ → 1⁺ → 2 → 2
@test Wcorr(
    theta1, phi1, theta2, phi2,
    State(0), M1(),
    State(1), Dipole(),
    State(2), Dipole(),
    State(2)) ≈ 0.0889651 rtol = 0.00001

# 0⁺ → 1⁻ → 2 → 2
@test Wcorr(
    theta1, phi1, theta2, phi2,
    State(0), E1(),
    State(1), Dipole(),
    State(2), Dipole(),
    State(2)) ≈ 0.077215 rtol = 0.00001

# 0⁺ → 2⁺ → 2 → 2
@test Wcorr(
    theta1, phi1, theta2, phi2,
    State(0), E2(),
    State(2), Dipole(),
    State(2), Dipole(),
    State(2)) ≈ 0.0562638 rtol = 0.00001
