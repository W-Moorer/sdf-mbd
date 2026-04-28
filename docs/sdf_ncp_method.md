# Signed Distance Field NCP Contact Method

## Motivation

Mesh-based contact can suffer from contact-pair switching, normal
discontinuity, and contact-search complexity. Penalty contact also introduces
stiffness sensitivity: the penetration error and stable time step depend on the
chosen penalty parameter. This prototype evaluates a classical, non-AI
alternative that combines signed distance fields with nonlinear
complementarity.

## SDF gap

For a contact point `x(q)` on a body and a signed distance field `phi`, the
normal gap is

```text
g(q) = phi(x(q)).
```

The convention in this prototype is `g > 0` outside the obstacle, `g = 0` on
the boundary, and `g < 0` in penetration. For the 2D ground plane `y = 0`,
`phi([x, y]) = y`.

## SDF normal

The contact normal is the normalized SDF gradient,

```text
n = grad phi(x(q)) / ||grad phi(x(q))||.
```

For exact analytical SDFs, `||grad phi|| = 1` away from nonsmooth points. The
implementation still normalizes defensively.

## SDF contact Jacobian

The SDF contact Jacobian maps generalized velocities to gap velocity:

```text
J_phi(q) = grad phi(x(q))^T * dx(q)/dq.
```

For the first point-mass prototype, `q = x`, so

```text
J_phi(q) = grad phi(q)^T.
```

The Python utility returns this row as a one-dimensional array of shape `(dim,)`.

## Signorini contact condition

Frictionless unilateral contact is modeled with the Signorini complementarity
condition,

```text
0 <= g(q),  lambda_n >= 0,  g(q) * lambda_n = 0.
```

Positive normal force is allowed only at closed contact, and open contact has
zero normal force.

## Smooth Fischer-Burmeister residual

The prototype enforces the complementarity equation with the smoothed
Fischer-Burmeister residual,

```text
Phi_eps(g, lambda) = sqrt(g^2 + lambda^2 + eps^2) - g - lambda.
```

The derivatives used by tests and future Newton assembly are

```text
dPhi/dg      = g / sqrt(g^2 + lambda^2 + eps^2) - 1
dPhi/dlambda = lambda / sqrt(g^2 + lambda^2 + eps^2) - 1.
```

As `eps -> 0`, this approaches the standard Fischer-Burmeister NCP function.
For `eps > 0`, the complementarity product is approximated rather than enforced
exactly.

## Point-mass implicit Euler residual

The first reproducible solver is a 2D point mass contacting an SDF plane. With

```text
q = [x, y]
v = [vx, vy]
M = m I
Q = [0, -m g0]
```

and unknown

```text
z = [vx_next, vy_next, lambda_n],
```

the step equations are

```text
q_next = q + dt * v_next
g_next = phi(q_next)
J_next = grad phi(q_next)^T
```

with residual

```text
R_v      = M (v_next - v) - dt * (Q + J_next^T lambda_n)
R_lambda = Phi_eps(g_next, lambda_n).
```

The implementation uses `scipy.optimize.root` when available and a small
finite-difference Newton fallback for this one-contact system.

## Penalty baseline

The baseline normal force is

```text
lambda_n = k_n * max(0, -g)
force    = lambda_n * n.
```

Optional normal damping can be added with `c_n * max(0, -g_dot)`. The first
paper scripts use the undamped baseline unless a script explicitly requests
damping.

## Current limitations

- The current prototype is frictionless.
- The current dynamics solver is limited to a 2D point mass and simple
  analytical SDF obstacles.
- Rigid-body rotational contact, multi-contact assembly, and patch contact are
  not yet included in the NCP solver.
- SDFs may be nonsmooth near medial axes, sharp corners, and min-composition
  boundaries.
- The smooth Fischer-Burmeister residual introduces an `eps`-dependent
  approximation.
- Deformable-body contact, frictional contact, and full Chrono C++ integration
  are future extensions.

## Future work

- Rigid-body SDF-NCP contact.
- Flexible multibody contact.
- Friction cone complementarity.
- Dynamic SDF reconstruction.
- Differentiable parameter identification.
- AI-assisted SDF correction may be studied later, but it is intentionally not
  part of this first classical paper prototype.

