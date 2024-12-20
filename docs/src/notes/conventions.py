import marimo

__generated_with = "0.10.6"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo

    import numpy as np
    import quaternion
    import sympy
    return mo, np, quaternion, sympy


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Three-dimensional space""")
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
def _(sympy, three_dimensional_coordinates):
    def three_dimensional_metric():
        (x, y, z), (r, θ, ϕ) = three_dimensional_coordinates()
        # The Cartesian metric is just the 3x3 identity matrix
        metric_cartesian = sympy.eye(3)
        # Compute the coordinate transformation to obtain the spherical metric
        jacobian = sympy.Matrix([x, y, z]).jacobian([r, θ, ϕ])
        metric_spherical = sympy.simplify((jacobian.T * metric_cartesian * jacobian))

        return metric_spherical

    three_dimensional_metric()
    return (three_dimensional_metric,)


@app.cell
def _(sympy, three_dimensional_coordinates):
    def three_dimensional_coordinate_basis():
        (x, y, z), (r, θ, ϕ) = three_dimensional_coordinates()
        basis = [
            [sympy.diff(xi, coord) for xi in (x, y, z)]
            for coord in (r, θ, ϕ)
        ]
        # Normalize each basis vector
        basis = [
            [sympy.simplify(comp / sympy.sqrt(sum(b**2 for b in vec))).subs(abs(sympy.sin(θ)), sympy.sin(θ))
            for comp in vec]
            for vec in basis
        ]

        return basis

    three_dimensional_coordinate_basis()
    return (three_dimensional_coordinate_basis,)


@app.cell
def _(sympy, three_dimensional_coordinates, three_dimensional_metric):
    def three_dimensional_volume_element():
        (x, y, z), (r, θ, ϕ) = three_dimensional_coordinates()
        metric_spherical = three_dimensional_metric()
        volume_element = sympy.simplify(sympy.sqrt(metric_spherical.det())).subs(abs(sympy.sin(θ)), sympy.sin(θ))

        return volume_element

    three_dimensional_volume_element()
    return (three_dimensional_volume_element,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Four-dimensional space""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### Geometric algebra""")
    return


@app.cell
def _():
    import galgebra

    from galgebra.ga import Ga
    return Ga, galgebra


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
def _(four_dimensional_coordinates, four_dimensional_volume_element):
    def Spin3_integral():
        import sympy
        (W, X, Y, Z), (R, α, β, γ) = four_dimensional_coordinates()
        volume_element = four_dimensional_volume_element()
        # Integrate over the unit sphere
        integral = (1/(16*sympy.pi**2)) * sympy.integrate(abs(sympy.sin(β)), (α, 0, 2*sympy.pi), (β, 0, 2*sympy.pi), (γ, 0, 2*sympy.pi))

        return integral

    Spin3_integral()
    return (Spin3_integral,)


@app.cell
def _(four_dimensional_coordinates, four_dimensional_volume_element):
    def SO3_integral():
        import sympy
        (W, X, Y, Z), (R, α, β, γ) = four_dimensional_coordinates()
        volume_element = four_dimensional_volume_element()
        # Integrate over the unit sphere
        integral = (1/(8*sympy.pi**2)) * sympy.integrate(sympy.sin(β), (α, 0, 2*sympy.pi), (β, 0, sympy.pi), (γ, 0, 2*sympy.pi))

        return integral

    SO3_integral()
    return (SO3_integral,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Rotations""")
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
