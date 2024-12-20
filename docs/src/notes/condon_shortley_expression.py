import marimo

__generated_with = "0.9.20"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def __(mo):
    mo.md(
        r"""
        Eq. (15) of Sec. 4³ (page 52) of Condon and Shortley (1935) defines the polar portion of the spherical harmonic function as

        \begin{equation}
        \Theta(\ell, m) = (-1)^\ell \sqrt{\frac{2\ell+1}{2} \frac{(\ell+m)!}{(\ell-m)!}}
        \frac{1}{2^\ell \ell!} \frac{1}{\sin^m\theta}
        \frac{d^{\ell-m}}{d(\cos\theta)^{\ell-m}} \sin^{2\ell}\theta.
        \end{equation}

        A footnote gives the first few values through $\ell=3$.  I explicitly test these explicit forms in [`SphericalFunctions.jl`](https://github.com/moble/SphericalFunctions.jl)`/test/conventions/condon_shortley.jl`.  Here, I want to verify that they are correct.

        Visually comparing, and accounting for some minor differences in simplification, I find that the expressions in the book are correct.  I also use the explicit expressions — as implemented in the test code and translated by AI — to check that sympy can simplify the difference to 0.
        """
    )
    return


@app.cell
def __():
    from IPython.display import display
    import marimo as mo
    import sympy
    from sympy import S

    θ = sympy.symbols("θ", real=True)

    def ϴ(ℓ, m):
        cosθ = sympy.symbols("cosθ", real=True)
        return (
            (-1)**ℓ
            * sympy.sqrt(
                ((2*ℓ+1) / 2)
                * (sympy.factorial(ℓ+m) / sympy.factorial(ℓ-m))
            )
            * (1 / (2**ℓ * sympy.factorial(ℓ)))
            * (1 / sympy.sin(θ)**m)
            #* sympy.diff(sympy.sin(θ)**(2*ℓ), sympy.cos(θ), ℓ-m)  # Can't differentiate wrt cos(θ), so we use a dummy and substitute
            * sympy.diff((1 - cosθ**2)**ℓ, cosθ, ℓ-m).subs(cosθ, sympy.cos(θ))
        ).simplify()
    return S, display, mo, sympy, Θ, θ


@app.cell
def __(S, display, Θ):
    for ℓ in range(4):
        for m in range(-ℓ, ℓ+1):
            display(ℓ, m, ϴ(S(ℓ), S(m)))
    return l, m


@app.cell
def __(S, sympy, Θ, θ):
    def compare_explicit_expression(ℓ, m):
        if (ℓ, m) == (0, 0):
            expression = sympy.sqrt(1/S(2))
        elif (ℓ, m) == (1, 0):
            expression = sympy.sqrt(3/S(2)) * sympy.cos(θ)
        elif (ℓ, m) == (2, 0):
            expression = sympy.sqrt(5/S(8)) * (2*sympy.cos(θ)**2 - sympy.sin(θ)**2)
        elif (ℓ, m) == (3, 0):
            expression = sympy.sqrt(7/S(8)) * (2*sympy.cos(θ)**3 - 3*sympy.cos(θ)*sympy.sin(θ)**2)
        elif (ℓ, m) == (1, 1):
            expression = -sympy.sqrt(3/S(4)) * sympy.sin(θ)
        elif (ℓ, m) == (1, -1):
            expression = sympy.sqrt(3/S(4)) * sympy.sin(θ)
        elif (ℓ, m) == (2, 1):
            expression = -sympy.sqrt(15/S(4)) * sympy.cos(θ) * sympy.sin(θ)
        elif (ℓ, m) == (2, -1):
            expression = sympy.sqrt(15/S(4)) * sympy.cos(θ) * sympy.sin(θ)
        elif (ℓ, m) == (3, 1):
            expression = -sympy.sqrt(21/S(32)) * (4*sympy.cos(θ)**2*sympy.sin(θ) - sympy.sin(θ)**3)
        elif (ℓ, m) == (3, -1):
            expression = sympy.sqrt(21/S(32)) * (4*sympy.cos(θ)**2*sympy.sin(θ) - sympy.sin(θ)**3)
        elif (ℓ, m) == (2, 2):
            expression = sympy.sqrt(15/S(16)) * sympy.sin(θ)**2
        elif (ℓ, m) == (2, -2):
            expression = sympy.sqrt(15/S(16)) * sympy.sin(θ)**2
        elif (ℓ, m) == (3, 2):
            expression = sympy.sqrt(105/S(16)) * sympy.cos(θ) * sympy.sin(θ)**2
        elif (ℓ, m) == (3, -2):
            expression = sympy.sqrt(105/S(16)) * sympy.cos(θ) * sympy.sin(θ)**2
        elif (ℓ, m) == (3, 3):
            expression = -sympy.sqrt(35/S(32)) * sympy.sin(θ)**3
        elif (ℓ, m) == (3, -3):
            expression = sympy.sqrt(35/S(32)) * sympy.sin(θ)**3
        else:
            raise ValueError(f"Unknown {ℓ=}, {m=}")
        return sympy.simplify(ϴ(S(ℓ), S(m)) - expression) == 0


    for _ℓ in range(4):
        for _m in range(-_ℓ, _ℓ+1):
            print((_ℓ, _m), "  \t", compare_explicit_expression(_ℓ, _m))

    return (compare_explicit_expression,)


if __name__ == "__main__":
    app.run()
