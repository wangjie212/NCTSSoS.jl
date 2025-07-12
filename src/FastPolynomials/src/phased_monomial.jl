struct PhasedMonomial
    phase::Tuple{Bool,Bool}
    mono::Monomial
end

function _merge_phase(p1::Tuple{Bool,Bool}, p2::Tuple{Bool,Bool})
    return (p1[1] ⊻ p2[1] ⊻ (p1[2] & p2[2]), p1[2] ⊻ p2[2])
end

function _phase_to_num(p1::Tuple{Bool,Bool})
    @match p1 begin
        (true, true) => -one(Float64) * im
        (true, false) => -one(ComplexF64)
        (false, true) => one(Float64) * im
        (false, false) => one(ComplexF64)
    end
end

function phased_monomial(m::Monomial)
    return PhasedMonomial((false, false), m)
end

# major difference at computing the constrain matrix
