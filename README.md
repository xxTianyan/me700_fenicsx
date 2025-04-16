# Linear Elastic Cantilever Beam Simulation (FEniCSx)

This project simulates the static deformation of a 3D cantilever beam under gravity using **linear elasticity theory** and finite element discretization via **FEniCSx**.

---

## 📐 Physical Problem Description

We consider a rectangular beam of length \( L \), width and height \( W \), fixed (clamped) at one end, and subject to a body force representing gravity:

- Domain: $\Omega = [0, L] \times [0, W] \times [0, W]$
- Boundary conditions:
  - **Clamped** at $x = 0 :  \mathbf{u} = \mathbf{0}$
  - **Free** elsewhere (no traction)
- Body force:

$$ \mathbf{f} = \begin{bmatrix} 0 \\\\ 0 \\\\ -\\rho g \end{bmatrix} $$

where $\rho$ is the density and $ g $ is a gravity-like scaling constant.

---

## 🧮 Governing Equations

We solve the **linear elasticity problem** in static equilibrium:

### 1. Strong form

$$
-\\nabla \cdot \boldsymbol{\sigma} = \mathbf{f} \quad \\text{in } \Omega
$$

with the **constitutive relation**:

$$
\boldsymbol{\sigma} = \\lambda \, (\nabla \cdot \mathbf{u}) \mathbf{I} + 2\mu \, \boldsymbol{\varepsilon}(\mathbf{u})
$$

and the **strain-displacement relation**:

$$
\boldsymbol{\varepsilon}(\mathbf{u}) = \frac{1}{2} \left( \nabla \mathbf{u} + (\nabla \mathbf{u})^T \\right)
$$

where:

- $\mathbf{u}$: displacement vector
- $\boldsymbol{\sigma}$: Cauchy stress tensor
- $\\lambda, \mu$: Lamé parameters (material constants)

---

## 🧩 Weak Formulation

The weak form of the equilibrium equation is:

$$
\\text{Find } \\mathbf{u} \\in V \text{ such that } \quad
a(\\mathbf{u}, \\mathbf{v}) = L(\\mathbf{v}) \quad \\forall \\mathbf{v} \\in V_0
$$

where

$$
a(\\mathbf{u}, \\mathbf{v}) = \int_\\Omega \boldsymbol{\sigma}(\mathbf{u}) : \boldsymbol{\varepsilon}(\mathbf{v}) \, d\\Omega
$$

$$
L(\\mathbf{v}) = \int_\\Omega \mathbf{f} \cdot \mathbf{v} \, d\\Omega
$$

---

## 🚗 Interpretation

The beam bends downward under its own weight. The maximum deflection is expected at the **free end**, and the displacement field is smooth and continuous due to the clamped constraint and linear elasticity assumption.

This setup is often used as a **benchmark problem** in solid mechanics and computational methods due to its simplicity and clear expected behavior.

## 📂 Files

- `beam_simulation.py` – Main script solving the elasticity problem
- `result.png` – Plot of z-displacement along the beam
