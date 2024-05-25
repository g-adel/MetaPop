function structVals(s)
    return Tuple(getfield(s, name) for name in fieldnames(typeof(s)))
end