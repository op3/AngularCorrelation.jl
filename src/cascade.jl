using EnumX
using HalfIntegers: Half, HalfInt

@enumx EMCharacter electric = 1 magnetic = -1 unknown = 0
@enumx Parity negative = -1 positive = 1 unknown = 0

"""
Representation of transition between two states
"""
struct Transition
    L::Int
    em_char::EMCharacter.T
    Lp::Int
    em_charp::EMCharacter.T
    δ::Real
end

"""
Invert the radiation character

char: original radiation character
"""
function alt_char(char::EMCharacter.T)
    if char ≡ Parity.negative
        Parity.positive
    elseif char ≡ Parity.positive
        Parity.negative
    else
        Parity.unknown
    end
end

E1(δ=0.0) = Transition(HalfInt(1), EMCharacter.electric, HalfInt(2), EMCharacter.magnetic, δ)
M1(δ=0.0) = Transition(HalfInt(1), EMCharacter.magnetic, HalfInt(2), EMCharacter.electric, δ)
E2(δ=0.0) = Transition(HalfInt(2), EMCharacter.electric, HalfInt(3), EMCharacter.magnetic, δ)
M2(δ=0.0) = Transition(HalfInt(2), EMCharacter.magnetic, HalfInt(3), EMCharacter.electric, δ)
Dipole(δ=0.0) = Transition(HalfInt(1), EMCharacter.unknown, HalfInt(2), EMCharacter.unknown, δ)
Quadrupole(δ=0.0) = Transition(HalfInt(2), EMCharacter.unknown, HalfInt(3), EMCharacter.unknown, δ)

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

@inline function check_em_char(parity_i::Half{T}, L::Int, em_char::EMCharacter.T, parity_f::Half{T}) where {T}
    if em_char ≡ EMCharacter.electric
        @assert parity_f ≡ parity_i^L
    elseif em_char ≡ EMCharacter.magnetic
        @assert parity_f ≡ parity_i^(L + 1)
    else
        true
    end
end

function check_transition(Si::State, γ::Transition, Sf::State)
    @assert Si.J + Sf.J ≥ γ.L ≥ abs(Si.J - Sf.J) "Selection rule violated for leading-order radiation"
    if γ.δ ≠ 0.0
        @assert Si.J + Sf.J ≥ γ.Lp ≥ abs(Si.J - Sf.J) "Selection rule violated for second-order radiation (γ.δ ≠ 0)"
    end
    @assert γ.L ≥ 1 "No monopole transitions"
    if Si.parity ≠ Parity.unknown && Sf.parity ≠ Parity.unknown
        if γ.em_char ≠ EMCharacter.unknown
            check_em_char(Si.parity, γ.L, γ.em_char, Sf.parity)
        end
        if γ.δ ≠ 0.0
            if γ.em_charp ≠ EMCharacter.unknown
                check_em_char(Si.parity, γ.Lp, γ.em_charp, Sf.parity)
            end
        end
    end
end

function check_cascade(cascade...)
    @assert length(cascade) % 2 ≡ 1
    @assert length(cascade) ≥ 3
    for i in 1:2:length(cascade)-4
        check_transition(cascade[i:i+2]...)
    end
end
