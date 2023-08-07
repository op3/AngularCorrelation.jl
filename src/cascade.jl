using EnumX
using HalfIntegers: Half, HalfInt

@enumx EMCharacter electric = 1 magnetic = -1 unknown = 0
@enumx Parity negative = -1 positive = 1 unknown = 0

"""
Radiation of a transition consisting of leading order and next-to-leading-order.

# Arguments:
- `L::Int`: Leading-order multipole
- `em_char::EMCharacter.T`: Leading-order radiation character
- `Lp::Int`: Next-to-leading-order multipole
- `em_charp::EMCharacter.T`: Next-to-leading-order radiation character
- `δ::Real`: Multipole mixing ratio between first to multipole orders (KSW convention)
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

# Arguments
- `char::EMCharacter.T`: original radiation character
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

"""
    E1(δ=0.0)

E1 radiation

# Arguments
- `δ`: Multipole mixing ration between E1 and M2 (KSW convention)
"""
E1(δ=0.0) = Transition(HalfInt(1), EMCharacter.electric, HalfInt(2), EMCharacter.magnetic, δ)

"""
    M1(δ=0.0)

M1 radiation

# Arguments
- `δ`: Multipole mixing ration between M1 and E2 (KSW convention)
"""
M1(δ=0.0) = Transition(HalfInt(1), EMCharacter.magnetic, HalfInt(2), EMCharacter.electric, δ)

"""
    E2(δ=0.0)

E2 radiation

# Arguments
- `δ`: Multipole mixing ration between E2 and M3 (KSW convention)
"""
E2(δ=0.0) = Transition(HalfInt(2), EMCharacter.electric, HalfInt(3), EMCharacter.magnetic, δ)

"""
    M2(δ=0.0)

M2 radiation

# Arguments
- `δ`: Multipole mixing ration between M2 and E3 (KSW convention)
"""
M2(δ=0.0) = Transition(HalfInt(2), EMCharacter.magnetic, HalfInt(3), EMCharacter.electric, δ)

"""
    Dipole(δ=0.0)

Dipole radiation with unknown radiation character

# Arguments
- `δ`: Multipole mixing ration between Dipole and Quadrupole (KSW convention)
"""
Dipole(δ=0.0) = Transition(HalfInt(1), EMCharacter.unknown, HalfInt(2), EMCharacter.unknown, δ)

"""
    Quadrupole(δ=0.0)

Quadrupole radiation with unknown radiation character

# Arguments
- `δ`: Multipole mixing ration between Quadrupole and Octupole (KSW convention)
"""
Quadrupole(δ=0.0) = Transition(HalfInt(2), EMCharacter.unknown, HalfInt(3), EMCharacter.unknown, δ)

Base.broadcastable(x::Transition) = Ref(x)

"""
A state represented by its quantum numbers
    
- `J::Half{T}`: Total angular momentum quantum number
- `parity::Parity.T`: Parity quantum number
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

struct CascadeException <: Exception
    msg::String
end

@inline function check_em_char(parity_i::Parity.T, L::Int, em_char::EMCharacter.T, parity_f::Parity.T)
    if em_char ≡ EMCharacter.electric
        if Int(parity_f) ≠ Int(parity_i) * (-1)^L
            throw(CascadeException("Wrong EMCharacter of radiation"))
        end
    elseif em_char ≡ EMCharacter.magnetic
        if Int(parity_f) ≠ Int(parity_i) * (-1)^(L + 1)
            throw(CascadeException("Wrong EMCharacter of radiation"))
        end
    else
        true
    end
end

function check_transition(Si::State, γ::Transition, Sf::State)
    if !(Si.J + Sf.J ≥ γ.L ≥ abs(Si.J - Sf.J))
        throw(CascadeException("Selection rule violated for leading-order radiation"))
    end
    if γ.δ ≠ 0.0 && !(Si.J + Sf.J ≥ γ.Lp ≥ abs(Si.J - Sf.J))
        throw(CascadeException("Selection rule violated for second-order radiation (γ.δ ≠ 0)"))
    end
    if γ.L ≡ 0
        throw(CascadeException("No monopole transitions"))
    end
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
    for i in 1:2:length(cascade)-2
        check_transition(cascade[i:i+2]...)
    end
end
