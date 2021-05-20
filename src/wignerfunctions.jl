
ϵ(m) = (m <= 0 ? 1 : (isodd(m) ? -1 : 1))


"""Return sign of input, with sign(0)=1"""
sign(m) = (m < 0 ? -1 : 1)


"""Compute Wigner's d matrix dˡₘₚ,ₘ(β)

Parameters
----------
expiβ : array_like
    Values of expi(i*β) on which to evaluate the d matrix.
out : array_like, optional
    Array into which the d values should be written.  It should be an array of
    floats, with size `self.dsize`.  If not present, the array will be created.
    In either case, the array will also be returned.
workspace : array_like, optional
    A working array like the one returned by Wigner.new_workspace().  If not
    present, this object's default workspace will be used.  Note that it is not
    safe to use the same workspace on multiple threads.

Returns
-------
d : array
    This is a 1-dimensional array of floats; see below.

See Also
--------
H : Compute a portion of the H matrix
D : Compute the full Wigner 𝔇 matrix
rotate : Avoid computing the full 𝔇 matrix and rotate modes directly
evaluate : Avoid computing the full 𝔇 matrix and evaluate modes directly

Notes
-----
This function is the preferred method of computing the d matrix for large ell
values.  In particular, above ell≈32 standard formulas become completely
unusable because of numerical instabilities and overflow.  This function uses
stable recursion methods instead, and should be usable beyond ell≈1000.

The result is returned in a 1-dimensional array ordered as

    [
        d(ell, mp, m, β)
        for ell in range(ell_max+1)
        for mp in range(-min(ℓ, mp_max), min(ℓ, mp_max)+1)
        for m in range(-ell, ell+1)
    ]

"""
function d!(d, w::WignerWorkspace, expiβ::Complex)
    ell_min = ℓₘᵢₙ(w.W)
    ell_max = ℓₘₐₓ(w.W)
    mp_max = m′ₘₐₓ(w.W)
    if mp_max < ell_max
        throw(DomainError("ℓₘₐₓ = $ℓₘₐₓ",
                          "Cannot compute full d matrix up to ℓₘₐₓ with m′ₘₐₓ only $(mp_max)"
        ))
    end

    H!(w, expiβ)

    for ℓ in ell_min:ell_max
        for m′ in -ℓ:ℓ
            for m in -ℓ:ℓ
                i_d = WignerDindex(ℓ, m′, m, ell_min)
                i_H = WignerHindex(ℓ, m′, m, mp_max)
                d[i_d] = ϵ(m′) * ϵ(-m) * w.Hwedge[i_H]
            end
        end
    end

    d
end


function d!(w::WignerWorkspace, β::Real)
    d = Array{T(w.W)}(undef, Wignerdsize(w.W))
    d!(w, expiβ, out)
end


"""
    D!(𝔇, w, R)
    D!(w, R)

Compute Wigner's 𝔇 matrix

Parameters
----------
𝔇 : array_like, optional
    Array into which the 𝔇 values should be written.  It should be an array of
    complex, with size `self.Dsize`.  If not present, the array will be
    created.  In either case, the array will also be returned.
workspace : optional
    A working array like the one returned by Wigner.new_workspace().  If not
    present, this object's default workspace will be used.  Note that it is not
    safe to use the same workspace on multiple threads.
R : Quaternion
    Array to be interpreted as a quaternionic array (thus its final dimension
    must have size 4), representing the rotations on which the 𝔇 matrix will be
    evaluated.

Returns
-------
D : array
    This is a 1-dimensional array of complex; see below.

See Also
--------
H : Compute a portion of the H matrix
d : Compute the full Wigner d matrix
rotate : Avoid computing the full 𝔇 matrix and rotate modes directly
evaluate : Avoid computing the full 𝔇 matrix and evaluate modes directly

Notes
-----
This function is the preferred method of computing the 𝔇 matrix for large ell
values.  In particular, above ell≈32 standard formulas become completely
unusable because of numerical instabilities and overflow.  This function uses
stable recursion methods instead, and should be usable beyond ell≈1000.

This function computes 𝔇ˡₘₚ,ₘ(R).  The result is returned in a 1-dimensional
array ordered as

    [
        𝔇(ell, mp, m, R)
        for ell in range(ell_max+1)
        for mp in range(-min(ℓ, mp_max), min(ℓ, mp_max)+1)
        for m in range(-ell, ell+1)
    ]

"""
function D!(𝔇, w::WignerWorkspace, R::Quaternion)
    ell_min = ℓₘᵢₙ(w.W)
    ell_max = ℓₘₐₓ(w.W)
    mp_max = m′ₘₐₓ(w.W)
    if mp_max < ell_max
        throw(DomainError("ℓₘₐₓ = $ℓₘₐₓ",
                          "Cannot compute full d matrix up to ℓₘₐₓ with m′ₘₐₓ only $(mp_max)"
        ))
    end

    to_euler_phases!(w.z, R)
    H!(w, w.z[2])
    complex_powers!(w.zₐpowers, w.z[1])
    complex_powers!(w.zᵧpowers, w.z[3])

    # 𝔇ˡₘₚ,ₘ(R) = dˡₘₚ,ₘ(R) exp[iϕₐ(m-mp)+iϕₛ(m+mp)] = dˡₘₚ,ₘ(R) exp[i(ϕₛ+ϕₐ)m+i(ϕₛ-ϕₐ)mp]
    # exp[iϕₛ] = R̂ₛ = hat(R[0] + 1j * R[3]) = zp
    # exp[iϕₐ] = R̂ₐ = hat(R[2] + 1j * R[1]) = zm.conjugate()
    # exp[i(ϕₛ+ϕₐ)] = zp * zm.conjugate() = z[2] = zᵧ
    # exp[i(ϕₛ-ϕₐ)] = zp * zm = z[0] = zₐ
    for ell in ell_min:ell_max
        for mp in -ell:-1
            i_D = WignerDindex(ell, mp, -ell, ell_min)
            for m in -ell:-1
                i_H = WignerHindex(ell, mp, m, mp_max)
                𝔇[i_D] = ϵ(mp) * ϵ(-m) * w.Hwedge[i_H] * conj(w.zᵧpowers[-m+1]) * conj(w.zₐpowers[-mp+1])
                i_D += 1
            end
            for m in 0:ell
                i_H = WignerHindex(ell, mp, m, mp_max)
                𝔇[i_D] = ϵ(mp) * ϵ(-m) * w.Hwedge[i_H] * w.zᵧpowers[m+1] * conj(w.zₐpowers[-mp+1])
                i_D += 1
            end
        end
        for mp in 0:ell
            i_D = WignerDindex(ell, mp, -ell, ell_min)
            for m in -ell:-1
                i_H = WignerHindex(ell, mp, m, mp_max)
                𝔇[i_D] = ϵ(mp) * ϵ(-m) * w.Hwedge[i_H] * conj(w.zᵧpowers[-m+1]) * w.zₐpowers[mp+1]
                i_D += 1
            end
            for m in 0:ell
                i_H = WignerHindex(ell, mp, m, mp_max)
                𝔇[i_D] = ϵ(mp) * ϵ(-m) * w.Hwedge[i_H] * w.zᵧpowers[m+1] * w.zₐpowers[mp+1]
                i_D += 1
            end
        end
    end
    𝔇
end


function D!(w::WignerWorkspace, R::Quaternion)
    𝔇 = Array{Complex{T(w.W)}}(undef, WignerDsize(w.W))
    D!(𝔇, w, R)
end
