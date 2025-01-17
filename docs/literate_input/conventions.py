import marimo

__generated_with = "0.10.6"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _():
    import marimo as mo

    import sympy
    # import numpy as np
    # import quaternion
    return mo, sympy


@app.cell(hide_code=True)
def _(mo):
    mo.md("""## Three-dimensional space""")
    return


@app.cell
def _():
    def three_dimensional_coordinates():
        import sympy
        # Make symbols representing spherical coordinates
        r, θ, ϕ = sympy.symbols("r θ ϕ", real=True, positive=True, zero=False)
        # Define Cartesian coordinates in terms of spherical
        x = r * sympy.sin(θ) * sympy.cos(ϕ)
        y = r * sympy.sin(θ) * sympy.sin(ϕ)
        z = r * sympy.cos(θ)
        return (x, y, z), (r, θ, ϕ)
    return (three_dimensional_coordinates,)


@app.cell
def _(mo, sympy, three_dimensional_coordinates):
    mo.output.append(mo.md("We define the spherical coordinates in terms of the Cartesian coordinates such that"))
    mo.output.append(
        [
            sympy.Eq(
                sympy.Symbol(f"{cartesian}"),
                spherical,
                evaluate=False
            )
            for cartesian, spherical in zip(("x", "y", "z"), three_dimensional_coordinates()[0])
        ]
    )
    return


@app.cell
def _(mo, three_dimensional_coordinates):
    def three_dimensional_metric():
        import sympy
        (x, y, z), (r, θ, ϕ) = three_dimensional_coordinates()
        # The Cartesian metric is just the 3x3 identity matrix
        metric_cartesian = sympy.eye(3)
        # Compute the coordinate transformation to obtain the spherical metric
        jacobian = sympy.Matrix([x, y, z]).jacobian([r, θ, ϕ])
        metric_spherical = sympy.simplify((jacobian.T * metric_cartesian * jacobian))

        return metric_spherical

    mo.output.append(mo.md("Using the Euclidean metric as the 3x3 identity in Cartesian coordinates, we can transform to spherical to obtain"))
    mo.output.append(three_dimensional_metric())
    return (three_dimensional_metric,)


@app.cell
def _(mo, three_dimensional_coordinates):
    def three_dimensional_coordinate_basis():
        import sympy
        (x, y, z), (r, θ, ϕ) = three_dimensional_coordinates()
        basis = [
            [sympy.diff(xi, coord) for xi in (x, y, z)]
            for coord in (r, θ, ϕ)
        ]
        # Normalize each basis vector
        basis = [
            sympy.Eq(
                sympy.Symbol(rf"\hat{{{symbol}}}"),
                sympy.Matrix([
                    sympy.simplify(comp / sympy.sqrt(sum(b**2 for b in vector))).subs(abs(sympy.sin(θ)), sympy.sin(θ))
                    for comp in vector
                ]),
                evaluate=False
            )
            for symbol, vector in zip((r, θ, ϕ), basis)
        ]

        return basis

    mo.output.append(mo.md(
        """
        Differentiating the $(x,y,z)$ coordinates with respect to $(r, θ, ϕ)$, we
        find the unit spherical coordinate basis vectors to have Cartesian components
        """
    ))
    mo.output.append(three_dimensional_coordinate_basis())
    return (three_dimensional_coordinate_basis,)


@app.cell
def _(mo, three_dimensional_coordinates, three_dimensional_metric):
    def three_dimensional_volume_element():
        import sympy
        (x, y, z), (r, θ, ϕ) = three_dimensional_coordinates()
        metric_spherical = three_dimensional_metric()
        volume_element = sympy.simplify(sympy.sqrt(metric_spherical.det())).subs(abs(sympy.sin(θ)), sympy.sin(θ))

        return volume_element


    mo.output.append(mo.md("The volume form in spherical coordinates is"))
    mo.output.append(three_dimensional_volume_element())
    return (three_dimensional_volume_element,)


@app.cell
def _(mo, three_dimensional_coordinates, three_dimensional_volume_element):
    def S2_normalized_volume_element():
        import sympy
        (x, y, z), (r, θ, ϕ) = three_dimensional_coordinates()
        volume_element = three_dimensional_volume_element().subs(r, 1)
        # Integrate over the unit sphere
        integral = sympy.integrate(sympy.sin(θ), (ϕ, 0, 2*sympy.pi), (θ, 0, sympy.pi))

        return volume_element / integral


    mo.output.append(mo.md(
        """
        Restricting to the unit sphere and normalizing so that the integral
        over the sphere is 1, we have the normalized volume element:
        """
    ))
    mo.output.append(S2_normalized_volume_element())
    return (S2_normalized_volume_element,)


@app.cell(hide_code=True)
def _(mo):
    mo.md("""## Four-dimensional space""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""### Geometric algebra""")
    return


@app.cell
def _():
    def three_dimensional_geometric_algebra():
        from galgebra.ga import Ga
        o3d, x, y, z = Ga.build("x y z", g=[1, 1, 1])

        i = z * y
        j = x * z
        k = y * x
        I = x * y * z

        return x, y, z, i, j, k, I
    return (three_dimensional_geometric_algebra,)


@app.cell
def _(three_dimensional_geometric_algebra):
    def check_basis_definitions():
        x, y, z, i, j, k, I = three_dimensional_geometric_algebra()
        return [
            i == z*y == -y*z,
            j == x*z == -z*x,
            k == y*x == -x*y,
            I == x*y*z == y*z*x == z*x*y == -x*z*y == -y*x*z == -z*y*x,
        ]

    check_basis_definitions()
    return (check_basis_definitions,)


@app.cell
def _(three_dimensional_geometric_algebra):
    def check_duals():
        x, y, z, i, j, k, I = three_dimensional_geometric_algebra()
        return [
            i == I.inv() * x,
            j == I.inv() * y,
            k == I.inv() * z,
        ]

    check_duals()
    return (check_duals,)


@app.cell(hide_code=True)
def _(mo):
    mo.md("""### Quaternions and Euler angles""")
    return


@app.cell
def _(three_dimensional_geometric_algebra):
    def check_basis_multiplication():
        x, y, z, i, j, k, I = three_dimensional_geometric_algebra()
        return [
            i*j == k,
            j*k == i,
            k*i == j,
            i*j*k == -1,
        ]

    check_basis_multiplication()
    return (check_basis_multiplication,)


@app.cell
def _():
    def four_dimensional_coordinates():
        import sympy
        # Make symbols representing spherical coordinates
        R, α, β, γ = sympy.symbols("R α β γ", real=True, positive=True)
        # Define Cartesian coordinates in terms of spherical
        W = R * sympy.cos(β/2) * sympy.cos((α + γ)/2)
        X = -R * sympy.sin(β/2) * sympy.sin((α - γ)/2)
        Y = R * sympy.sin(β/2) * sympy.cos((α - γ)/2)
        Z = R * sympy.cos(β/2) * sympy.sin((α + γ)/2)

        return (W, X, Y, Z), (R, α, β, γ)
    return (four_dimensional_coordinates,)


@app.cell
def _(four_dimensional_coordinates):
    def four_dimensional_metric():
        import sympy
        (W, X, Y, Z), (R, α, β, γ) = four_dimensional_coordinates()
        # The Cartesian metric is just the 4x4 identity matrix
        metric_cartesian = sympy.eye(4)
        # Compute the coordinate transformation to obtain the spherical metric
        jacobian = sympy.Matrix([W, X, Y, Z]).jacobian([R, α, β, γ])
        metric_spherical = sympy.simplify((jacobian.T * metric_cartesian * jacobian))

        return metric_spherical

    four_dimensional_metric()
    return (four_dimensional_metric,)


@app.cell
def _(four_dimensional_metric):
    def four_dimensional_volume_element():
        import sympy
        metric_spherical = four_dimensional_metric()
        # Compute the volume element
        volume_element = sympy.sqrt(metric_spherical.det()).simplify()

        return volume_element

    four_dimensional_volume_element()
    return (four_dimensional_volume_element,)


@app.cell
def _(four_dimensional_coordinates, four_dimensional_volume_element, mo):
    def Spin3_normalized_volume_element():
        import sympy
        (W, X, Y, Z), (R, α, β, γ) = four_dimensional_coordinates()
        volume_element = four_dimensional_volume_element().subs(R, 1)
        # Integrate over the unit sphere
        integral = sympy.integrate(volume_element, (α, 0, 2*sympy.pi), (β, 0, 2*sympy.pi), (γ, 0, 2*sympy.pi))

        return volume_element / integral


    mo.output.append(mo.md(
        r"""
        Restricting to the unit three-sphere, $\mathrm{Spin}(3)$, and normalizing so that the integral
        over the sphere is 1, we have the normalized volume element:
        """
    ))
    mo.output.append(Spin3_normalized_volume_element())
    return (Spin3_normalized_volume_element,)


@app.cell
def _(four_dimensional_coordinates, four_dimensional_volume_element, mo):
    def SO3_normalized_volume_element():
        import sympy
        (W, X, Y, Z), (R, α, β, γ) = four_dimensional_coordinates()
        volume_element = four_dimensional_volume_element().subs(R, 1)
        # Integrate over the unit sphere
        integral = sympy.integrate(volume_element, (α, 0, 2*sympy.pi), (β, 0, sympy.pi), (γ, 0, 2*sympy.pi))

        return volume_element / integral


    mo.output.append(mo.md(
        r"""
        Restricting to $\mathrm{SO}(3)$ and normalizing so that the integral
        over the sphere is 1, we have the normalized volume element:
        """
    ))
    mo.output.append(SO3_normalized_volume_element())
    return (SO3_normalized_volume_element,)


@app.cell(hide_code=True)
def _(mo):
    mo.md("""## Rotations""")
    return


@app.cell
def _():
    # Check that Euler angles (α, β, γ) rotate the basis vectors onto (r, \theta, \phi)


    return


if __name__ == "__main__":
    app.run()
