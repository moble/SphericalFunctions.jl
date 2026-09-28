# Serialization writes every pointer as a null pointer.  An FFTW plan object wraps a pointer
# to a plan in the memory of the process that made it, so the plans of a transform would
# arrive in another process — or be read back from a file — as plans that crash the process
# when they are executed.  The plans of a `RingPlans` are therefore not reconstructed from
# what was written: they are made again in the receiving process, of the sizes and with the
# planner options that were written.  What is written is the default serialization of the
# struct, its fields in order, and the plans it contains are read and discarded.  (A plan of
# the generic FFT holds no pointer, and would survive, but is remade in the same way.)
function Serialization.deserialize(
    s::Serialization.AbstractSerializer, ::Type{RingPlans{T, P, BP}}
) where {T, P, BP}
    sizes = Serialization.deserialize(s)
    index = Serialization.deserialize(s)
    flags = Serialization.deserialize(s)
    timelimit = Serialization.deserialize(s)
    Serialization.deserialize(s)  # the forward plans of the sending process
    Serialization.deserialize(s)  # and the backward plans
    ring_plans(T, sizes, index, flags, timelimit)::RingPlans{T, P, BP}
end
