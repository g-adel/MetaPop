mutable struct SA
    a
end

mutable struct SB
    sa::SA
end

sa = SA(1)
sb = SB(sa)
sa.a=2
@show sa, sb
