using AngularCorrelation
using Test

@test E1() == Transition(1, EMCharacter.electric, 2, EMCharacter.magnetic, 0.0)
@test M1() == Transition(1, EMCharacter.magnetic, 2, EMCharacter.electric, 0.0)
@test Dipole() == Transition(1, EMCharacter.unknown, 2, EMCharacter.unknown, 0.0)
@test E2() == Transition(2, EMCharacter.electric, 3, EMCharacter.magnetic, 0.0)
@test M2() == Transition(2, EMCharacter.magnetic, 3, EMCharacter.electric, 0.0)
@test Quadrupole() == Transition(2, EMCharacter.unknown, 3, EMCharacter.unknown, 0.0)
@test E1(3.3) == Transition(1, EMCharacter.electric, 2, EMCharacter.magnetic, 3.3)
@test M1(-3.3) == Transition(1, EMCharacter.magnetic, 2, EMCharacter.electric, -3.3)
@test Dipole(2.2) == Transition(1, EMCharacter.unknown, 2, EMCharacter.unknown, 2.2)
@test E2(-1.0) == Transition(2, EMCharacter.electric, 3, EMCharacter.magnetic, -1.0)
@test M2(1.0) == Transition(2, EMCharacter.magnetic, 3, EMCharacter.electric, 1.0)
@test Quadrupole(2.2) == Transition(2, EMCharacter.unknown, 3, EMCharacter.unknown, 2.2)
