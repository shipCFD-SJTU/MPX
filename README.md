# MPX — Multi-Physics eXtended Solver

**An OpenFOAM-based research solver for temperature–salinity-driven density-stratified flows.**

MPX (*Multi-Physics eXtended Solver*) extends the OpenFOAM framework to study fluid motion in which **temperature (T)** and **salinity (S)** jointly determine density, buoyancy, and stratification. The project is motivated by ocean and marine hydrodynamics, including the interaction of stratified environments with submerged bodies, wakes, internal waves, and moving structures.

MPX is developed as a **separately maintained OpenFOAM extension**: it aims to minimize modifications to the upstream OpenFOAM source tree while retaining access to its finite-volume infrastructure and numerical capabilities.

> **Project status:** Research software under active development. The initial public code is based on **OpenFOAM v2606**. The supported solvers, initialization options, and dynamic-mesh capabilities depend on the contents of the checked-out revision. Please consult the source code and example cases before using a feature in production.

## Contents

- [Why MPX?](#why-mpx)
- [Key capabilities](#key-capabilities)
- [Governing equations](#governing-equations)
- [Requirements](#requirements)
- [Installation and compilation](#installation-and-compilation)
- [Running a case](#running-a-case)
- [Verification and validation](#verification-and-validation)
- [Project structure](#project-structure)
- [Roadmap](#roadmap)
- [Contributing](#contributing)
- [Citation](#citation)
- [License and third-party software](#license-and-third-party-software)
- [Acknowledgments](#acknowledgments)

## Why MPX?

In the ocean, density variations are commonly controlled by both temperature and salinity. Representing density using only one transported scalar can be insufficient when the objective is to investigate thermohaline transport, mixing, or the evolution of a stratified wake.

MPX therefore treats temperature and salinity as **separate transported fields** and evaluates density using a configurable equation of state (EOS). This makes it possible to investigate an important physical question:

> **Identical initial density stratification does not necessarily imply identical temperature–salinity transport, mixing, or wake evolution.**

The framework is intended to support both numerical-method development and fundamental investigations of stratified-flow physics.

## Key capabilities

| Area                        | Capability                                                   |
| --------------------------- | ------------------------------------------------------------ |
| CFD framework               | OpenFOAM finite-volume discretization and pressure–velocity coupling |
| Stratification              | Temperature–salinity-driven density variations               |
| Scalar transport            | Separate transport equations for temperature and salinity    |
| Buoyancy                    | Boussinesq approximation for variable-density buoyancy effects |
| Initial profiles            | Linear and hyperbolic-tangent (`tanh`) stratification; configuration of profiles from measured T/S data |
| Density modeling            | Linear temperature–salinity EOS and TEOS-10-based seawater-density modeling |
| Marine applications         | Stratified flow past spheres and underwater bodies; wakes, internal waves, and mixing |
| Motion and mesh integration | Development and integration pathways for six-degree-of-freedom (6DoF) motion and overset meshes |
| Code organization           | External development and compilation relative to the main OpenFOAM source tree |

**Important distinctions:**

- A *linear* or *tanh* profile describes the **initial stratification**, not the equation of state itself.
- The actual 6DoF/overset workflow must be checked against the selected solver and case; these capabilities should not be assumed available in every build.
- For a TEOS-10-based EOS, the temperature, salinity, and pressure conventions and units must be documented and used consistently.

## Governing equations

MPX is developed around the incompressible/Boussinesq formulation. In schematic form, the governing equations are:

**Mass conservation**

$$
\nabla \cdot \mathbf{U} = 0.
$$

**Momentum conservation**

$$
\rho_0\left(\frac{\partial \mathbf{U}}{\partial t}
+ \mathbf{U}\cdot\nabla\mathbf{U}\right)
= -\nabla p + \nabla\cdot\boldsymbol{\tau}
+ (\rho-\rho_0)\mathbf{g},
$$

where the hydrostatic reference state is absorbed into the pressure definition.

**Temperature transport**

$$
\frac{\partial T}{\partial t}
+ \nabla\cdot(\mathbf{U}T)
= \nabla\cdot(\alpha_{T,\mathrm{eff}}\nabla T).
$$

**Salinity transport**

$$
\frac{\partial S}{\partial t}
+ \nabla\cdot(\mathbf{U}S)
= \nabla\cdot(\alpha_{S,\mathrm{eff}}\nabla S).
$$

**Density closure**

$$
\rho = \rho(T,S,p) \quad \text{(according to the selected EOS)}.
$$

A commonly used linearized EOS is

$$
\rho = \rho_0\left[1-\beta_T(T-T_0)+\beta_S(S-S_0)\right],
$$

where $\beta_T$ and $\beta_S$ are the thermal expansion and haline contraction coefficients. The reference density $\rho_0$ is used in the Boussinesq inertia terms; density differences enter the buoyancy contribution.

The equations above describe the physical model, **not a line-by-line specification of every discretization or turbulence closure** in the repository. For the implementation details, consult the solver source, EOS routines, and OpenFOAM case dictionaries.

## Requirements

The initial development environment is:

| Component                | Reference environment                                        |
| ------------------------ | ------------------------------------------------------------ |
| Operating system         | Ubuntu 24.04 LTS (native Linux or WSL2)                      |
| CFD platform             | **OpenFOAM v2606** (the [openfoam.com](https://www.openfoam.com/) distribution) |
| Compiler and build tools | A working compiler toolchain and OpenFOAM `wmake`            |
| Version control          | Git                                                          |
| Visualization (optional) | ParaView                                                     |

Other OpenFOAM releases have **not** been established as compatible by this README. A successful compilation on v2606 does not guarantee compatibility with releases using different APIs.

Install and configure OpenFOAM first. Official documentation: [OpenFOAM v2606](https://doc.openfoam.com/2606/) and [source packages](https://dl.openfoam.com/source/v2606/).

## Installation and compilation

### 1. Load the OpenFOAM environment

Adjust the path to your local installation:

```bash
source "$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc"
```

Check the active environment:

```bash
echo "$WM_PROJECT_VERSION"
echo "$WM_PROJECT_DIR"
command -v wmake
```

The active OpenFOAM version should be `v2606` (or `2606`, depending on the environment).

### 2. Obtain the MPX source code

Replace the following GitHub placeholders with this repository's actual owner and repository name:

```bash
git clone https://github.com/shipCFD-SJTU/MPX
cd MPX
```

### 3. Compile MPX

MPX is intended to be compiled **outside** the upstream OpenFOAM source tree.

If the repository contains a top-level `Allwmake` script, build with:

```bash
chmod +x Allwmake
./Allwmake
```

If there is **no** top-level `Allwmake`, enter the relevant solver or utility source directory containing `Make/files` and `Make/options`, and compile that component with:

```bash
wmake
```

The correct source directory and build order depend on the files included in this repository. Do not run `wmake` blindly in the repository root if no corresponding build configuration is present.

### 4. Verify the build

Check the target executable names in each component's `Make/files`. If the respective targets are named `MPX` and `initMPXStratification`, for example, verify them with:

```bash
command -v MPX
command -v initMPXStratification
```

If an executable is not found, check the compiler output, `Make/files` target path, and OpenFOAM user application directory (`$FOAM_USER_APPBIN`).

## Running a case

An MPX case follows the usual OpenFOAM organization, typically including:

```text
case/
├── 0/                  # Initial and boundary fields (e.g., U, T, S, pressure)
├── constant/           # Physical properties, mesh, and model dictionaries
└── system/             # Solver settings and discretization controls
```

The exact field names and dictionaries are determined by the solver revision and supplied example cases.

A typical workflow is:

1. Generate or import the computational mesh.
2. Configure the temperature and salinity fields, their transport properties, and the selected EOS.
3. Initialize the prescribed T/S stratification, checking that the density profile is physically stable.
4. Select boundary conditions, numerical schemes, time-step controls, and any turbulence model.
5. Run the compiled MPX solver and inspect residuals, mass conservation, and changes in the background stratification.
6. Post-process velocity, temperature, salinity, density, wake structure, and other relevant flow quantities.

**Illustrative commands** — use only when the corresponding files and executables exist:

```bash
cd /path/to/your/case

# Only for cases configured with system/blockMeshDict:
blockMesh

# Only if the initializer is built and configured for this case:
initMPXStratification

# Only if the compiled solver executable is named MPX:
MPX
```

For reproducibility, record the Git commit, OpenFOAM build, EOS parameters, transport coefficients, mesh, initial T/S profiles, and numerical settings with each simulation.

## Verification and validation

The MPX verification and validation program emphasizes both **numerical accuracy** and **the preservation of an intended density-stratified environment**.

### Background stratification

- Verify temperature, salinity, and density profiles against their prescribed initial distributions.
- Check density drift, spurious currents, and long-duration stability in nominally quiescent stratified fluid.
- Assess sensitivity to mesh resolution, time step, numerical schemes, and scalar diffusivities.
- Test linear, `tanh`, and data-prescribed profiles separately from the choice of EOS.

### Flow past a sphere

A reference test problem in the MPX research program is a sphere in linearly stratified flow:

| Parameter                    | Reference value |
| ---------------------------- | --------------- |
| Sphere diameter, $D$         | **0.5 m**       |
| Reynolds number, $Re=UD/\nu$ | **3700**        |
| Froude number, $Fr=U/(ND)$   | **1.0**         |
| Stratification               | Linear          |

Here, $U$ is the reference inflow speed, $\nu$ is kinematic viscosity, and $N$ is buoyancy frequency, with $N^2=-(g/\rho_0)\,d\rho/dz$ when the vertical coordinate $z$ points upward.

This table defines a **reference configuration**, not a claim of a particular measured error or a completed validation dataset. Quantitative comparisons, numerical uncertainty, and experimental references should be reported in accompanying case documentation or publications.

Further applications of interest include more complex underwater bodies, thermohaline wakes, internal-wave generation, and coupled body motion.

## Project structure

MPX is organized around the separation of the OpenFOAM platform from MPX-specific numerical methods. Depending on the repository revision, components may include:

- **Solver code:** pressure–velocity coupling and integration of the T/S transport equations.
- **Stratification initialization utilities:** creation of specified initial temperature/salinity distributions.
- **EOS components:** evaluation of density from the thermodynamic state.
- **Example cases and validation data:** reproducible configurations for verification and research.
- **Documentation:** build instructions, numerical methodology, and usage notes.

## Roadmap

Areas targeted for continued development include:

- [ ] Broaden automated verification and regression tests.
- [ ] Publish reproducible stratification-preservation and sphere benchmarks.
- [ ] Improve configuration, documentation, and testing of the available EOS options.
- [ ] Extend and validate 6DoF and overset-mesh coupling workflows.
- [ ] Test additional OpenFOAM versions without modifying their core source trees.
- [ ] Add more canonical cases for underwater-body wakes and ocean stratification.

The checklist describes **planned repository milestones**; it does not downgrade or supersede research capabilities already explored in development branches.

## Contributing

Scientific discussion, reproducible bug reports, documentation improvements, and contributions are welcome through GitHub Issues and Pull Requests.

When reporting a problem, please include:

- The MPX commit or release tag and the OpenFOAM version.
- Operating system, compiler, and build options.
- A minimal case or the relevant configuration dictionaries.
- Build output or solver logs and the steps needed to reproduce the issue.

Please avoid submitting confidential research data or third-party code without permission.

## Citation

If MPX contributes to your research, please acknowledge the software and identify the **exact repository revision or released version** used in your simulations. A formal `CITATION.cff`, archival DOI, and method-paper citation can be added when those resources are available.

A suggested interim acknowledgment is:

> The simulations were performed using MPX (Multi-Physics eXtended Solver), an OpenFOAM-based research solver for temperature–salinity-driven density-stratified flows. https://github.com/shipCFD-SJTU/MPX

## License and third-party software

MPX builds on **OpenFOAM**, which is distributed under the GNU General Public License. See the [OpenFOAM licensing information](https://www.openfoam.com/licensing) and the license notices in the upstream source code.

**For the MPX code itself, consult this repository's `LICENSE` file and individual source-file headers.** If no license has been published, the availability of source code on GitHub alone does not establish permission to copy, modify, or redistribute it. Any OpenFOAM-derived code and third-party dependencies must retain their applicable notices and comply with their respective licenses.

MPX is an independent research project; it is not presented as an official OpenFOAM release or as being endorsed by the OpenFOAM developers.

## Acknowledgments

MPX benefits from the open-source OpenFOAM ecosystem and from the broader research community working on computational fluid dynamics, buoyancy-driven flows, and ocean stratification.

---

**MPX — Multi-Physics eXtended Solver**  
*Temperature–salinity physics · Density stratification · Marine hydrodynamics*
