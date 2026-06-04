# Candidate for upstreaming to ThreeBodyDecays (coupling_schemes.jl).

"""
    possible_ls_more(jp1, jp2; jp)

Allowed `(two_l, two_s)` for a two-body vertex. Same spin ranges as
`ThreeBodyDecays.possible_ls`. Parity is enforced with the usual `⊗` rule unless
any line has [`UndefinedParity`](@ref) (`'±'`).
"""
function possible_ls_more(jp1::SpinParity, jp2::SpinParity; jp::SpinParity)
    two_ls = Vector{Tuple{Int,Int}}(undef, 0)
    for two_s in abs(jp1.two_j - jp2.two_j):2:(jp1.two_j+jp2.two_j)
        for two_l in abs(jp.two_j - two_s):2:(jp.two_j+two_s)
            if jp1.p == UndefinedParity ||
               jp2.p == UndefinedParity ||
               jp.p == UndefinedParity
                push!(two_ls, (two_l, two_s))
            elseif jp1.p ⊗ jp2.p ⊗ jp.p == (isodd(div(two_l, 2)) ? '-' : '+')
                push!(two_ls, (two_l, two_s))
            end
        end
    end
    return sort(two_ls, by=x -> x[1])
end
