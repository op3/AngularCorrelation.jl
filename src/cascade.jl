using EnumX
using HalfIntegers: Half, HalfInt

@enumx EMCharacter electric = -1 magnetic = 1 unknown = 0
@enumx Parity negative = -1 positive = 1 unknown = 0

"""
Representation of transition between two states
"""
struct Transition
    L::UInt
    em_char::EMCharacter.T
    Lp::UInt
    em_charp::EMCharacter.T
    delta::Real
end

"""
Invert the radiation character

char: original radiation character
"""
function alt_char(char::EMCharacter.T)
    if char == Parity.negative
        Parity.positive
    elseif char == Parity.positive
        Parity.negative
    else
        Parity.unknown
    end
end

E1(delta=0.0) = Transition(HalfInt(1), EMCharacter.electric, HalfInt(2), EMCharacter.magnetic, delta)
M1(delta=0.0) = Transition(HalfInt(1), EMCharacter.magnetic, HalfInt(2), EMCharacter.electric, delta)
E2(delta=0.0) = Transition(HalfInt(2), EMCharacter.electric, HalfInt(3), EMCharacter.magnetic, delta)
M2(delta=0.0) = Transition(HalfInt(2), EMCharacter.magnetic, HalfInt(3), EMCharacter.electric, delta)
Dipole(delta=0.0) = Transition(HalfInt(1), EMCharacter.unknown, HalfInt(2), EMCharacter.unknown, delta)
Quadrupole(delta=0.0) = Transition(HalfInt(2), EMCharacter.unknown, HalfInt(3), EMCharacter.unknown, delta)

Base.broadcastable(x::Transition) = Ref(x)

"""
Representation of quantum numbers of a state
"""
struct State{T}
    J::Half{T}
    parity::Parity.T
end

Base.broadcastable(x::State) = Ref(x)

function State(J::Int, parity::Parity.T=Parity.unknown)
    State(HalfInt(J), parity)
end

function State(J::Rational, parity::Parity.T=Parity.unknown)
    State(HalfInt(J), parity)
end
