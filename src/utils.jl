function tight_vector_pyconvert(pyv)
    return pyconvert(String, pyv.dtype.name) == "object" ? string.(collect(pyv)) : pyconvert(Vector, pyv)
end

to_dataframe(pd) = DataFrame([pyconvert(String, col) => tight_vector_pyconvert(pd[col]) for col in pd.columns])
