struct PhasedMonomial
    phase::Tuple{Bool,Bool}
    mono::Monomial
end

# major difference at computing the constrain matrix
