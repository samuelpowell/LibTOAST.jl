# libTOAST.jl


[![Build Status](https://travis-ci.org/samuelpowell/libTOAST.jl.svg?branch=master)](https://travis-ci.org/samuelpowell/libTOAST.jl)

[![Coverage Status](https://coveralls.io/repos/samuelpowell/libTOAST.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/samuelpowell/libTOAST.jl?branch=master)

[![codecov.io](http://codecov.io/github/samuelpowell/libTOAST.jl/coverage.svg?branch=master)](http://codecov.io/github/samuelpowell/libTOAST.jl?branch=master)


libTOAST.jl is a low-level interface to the [TOAST++](https://github.com/toastpp/toastpp) library. Whilst TOAST++ is designed as an end-to-end solution for forward modelling and image reconstruction in Diffuse Optical Tomography (DOT), the libTOAST.jl interface provides an (opinionated) interface to a subset of its underlying finite-element, raster mapping, and regularisation functionality. DOT specific functionality is implemented in a separate Julia package.

libTOAST.jl allows one to easily:

1. solve certain types of second order partial differential equations in two- and three-dimensions using the (Galerkin) Finite Element Method;
2. map functions between unstructured meshes and alternative bases (e.g. pixels, voxels) defined in a uniform rasterisation;
3. evaluate the value and derivatives of a number of spatial regularisation functionals.

## Compatibilty

libTOAST.jl directly calls the underlying C++ library using [Cxx.jl](https://github.com/Keno/Cxx.jl) which is currently supported only on Linux and macOS. Windows support is *expected* with the release of Julia 0.6.

## Installation

Ensure that you have a working Cxx.jl installation, then install and build libTOAST.jl package as follows. Installation on Linux requires that you have the `patchelf` utility installed.

```
julia> ] add git@github.com:samuelpowell/libTOAST.jl.git
julia> ] build libTOAST


The libTOAST.jl build process will download the requisite binaries and the header files required to call the library, these are stored in the package directory. You can test that libTOAST.jl is properly configured by running the unit tests.

```
Pkg.test("libTOAST")
```

## A first example

By means of an introduction, we will use libTOAST.jl to solve a lossy diffusion equation over a domain Ω

```
[-∇⋅κ(r)∇ + μ(r)] ϕ(r) = δ(r) (r ∈ Ω),
```

with a Robin condition on the boundary δΩ

```
ϕ(r) + γ n⋅∇ϕ(r) = 0 (r ∈ δΩ),
```

where `n` is the unit outward normal to the boundary.

Solving this PDE numerically using the Finite Element Method can be achieved as follows:

1. find the weak form of the equation through multiplication by a test function and integration by parts;
2. subdivide the domain into a mesh of non-overlapping elements joined at N vertex nodes, and define a set of basis functions over this domain;
3. choose the test functions in the weak formulation to be the same as the aformentioned basis in the mesh;
4. approximate the solution and the parameters parameters in the same basis, e.g,

```
ϕ(r) ≈ ϕʰ(r) = ∑ᵢ ϕᵢ uᵢ(r).
```

We further assume a representation of the parameters γ(r), κ(r), and μ(r) in the same basis. Performing this process results in a linear system of equations

```
S ϕ = Q,
```

where the system matrix `S = K + M + F` is calculated as

```
Kᵢⱼ = ∑ᵥ κᵥ ∫_Ω  uᵥ(r) ∇uᵢ(r)⋅∇uⱼ(r) dr,
Mᵢⱼ = ∑ᵥ μᵥ ∫_Ω  uᵥ(r) uᵢ(r) uⱼ(r) dr,
Fᵢⱼ = ∑ᵥ γᵥ ∫_δΩ uᵥ(r) uᵢ(r) uⱼ(r) dr,
```
the source term

```
Qᵢ = ∑ᵢ qᵢ ∫_Ω uᵢ(r) dr,
```

and `ϕ` is the vector of nodal coefficients that specified the solution in the chosen basis.

Items (1) - (4) are typically performed by hand. The core utility of the finite-element functionality of libTOAST.jl is to manage the mesh, compute the integrals of various products of basis functions, and assemble the matrices required to find the solution vector. We will demonstrate this process in the following five steps. The full source code for this example is located in `/examples/readme.jl`.

### 1. Load the TOAST module

To begin, load the TOAST module.

```
julia> using TOAST
Default thread count set to 8
```

The TOAST++ backend reports that it has initialised a pool of eight threads in which it will perform tasks such as system matrix assembly.

### 2. Build the mesh

In this example we will define the domain by loading a mesh from a file.

```
julia> meshdir = joinpath(libTOAST._jl_toast_dir, "test", "2D", "meshes")

julia> mesh = Mesh(joinpath(meshdir, "circle25_32.msh"));
INFO: Loaded 250114 bytes from /Users/spowell/.julia/v0.5/TOAST/deps/toastpp-2.0.2/test/2D/meshes/circle25_32.msh.
INFO: Initiliasing mesh.
INFO: Calculating sparsity pattern.
INFO: Number of vertices: 3511, number of nonzeros: 24211
```

This mesh defines a circular two-dimensional domain of radius 25mm.

libTOAST.jl reports that it has loaded the mesh. In TOAST++ the mesh and the finite-element subspace are tightly coupled: a description of the mesh includes a specification of the finite-element shape functions; as such libTOAST.jl is able to precompute the sparsity pattern which will be used for system matrix assembly.

### 3. Define the parameters

There are three parameters in the problem, κ(r), μ(r), and γ(r). We directly provide the discrete form of the parameters as a set of nodal coefficients in the same basis as the shape functions used for the finite-element solution.

Functions in nodal basis are represented in `NodalCoeff` types, which hold a reference to the `mesh`.

```
julia> μ = fill(NodalCoeff, mesh, 0.01);
julia> κ = fill(NodalCoeff, mesh, 0.3);
julia> γ = fill(NodalCoeff, mesh, 0.5);
```

### 4. Build the system matrix

To build the system matrix we first instantiate a `SystemMatrix`

```
julia> S = SystemMatrix(mesh)
```

which is a thin wrapper around a compressed sparse row (CSR) matrix used by TOAST. Since Julia's native sparse arrays are compressed sparse column (CSC) libTOAST.jl, uses this wrapper until the final matrix is ready for use, this allows us to add components of the system matrix in-place without additional allocations.

We now assemble and add each of the components to the system matrix.

```
julia> assemble!(S, PDD, κ);
julia> assemble!(S, PFF, μ);
julia> assemble!(S, BNDPFF, γ);
```

The second input to the `assemble!` method defines the integral to be performed. The available integrals are enumerated later, but for now it will suffice to understand that `P` represents a parameter in the mesh basis, `F` represents a basis function, `D` represents a derivative of a basis function, and the `BND` prefix indicates the integration is only over the boundary.

When we are finally ready to use the system matrix in Julia, we can do so as follows.

```
sysmat = sparse(S)
```

### 4. Build a source term

Equipped with a system matrix describing the discrete form of the PDE we can solve the system against an arbitrary source term (right hand side). In general, one may build arbitrary source terms using the general-purpose assembly routines. In this case, however, we choose to solve against a point source located at the centre of the domain. For this purposes libTOAST.jl exposes a specific type and special-purpose method `assemble` method,

```
julia> source = PointSource(mesh, [0.,0.], Isotropic)
julia> Q = assemble(source)
```

which builds a source vector representing an isotropic unitary point source centred at the specified co-ordinated.

### 5. Solve the system (and visualise)

The resultant system of equations can be solved using your method of preference. For example, we can simply rely upon Julia's backslash operator to determine the properties of the system matrix and undertake a suitable direct solve.

```
julia> ϕ = sparse(S)\Q
```

If you have a working installation of PyPlot, you can easily visualise the results by extracting the vertices from the mesh and plotting the triangulation.

```
julia> using PyPlot
julia> vtx, = data(mesh);
julia> plot_trisurf(vtx[:,1], vtx[:,2], ϕ)
```

## Concepts, Types, and Methods

### Meshes and the FE subspace

In TOAST++ the mesh and the underlying finite-element subspace are tightly coupled: the type of the finite-element is chosen during initialisation of the mesh. libTOAST.jl only supports homogenous meshes of piecewise-linear triangular (2D) or tetrahedral (3D) meshes. Element-wise assembly is a planned feature.

Meshes can be specified by supplying two parameters:

1. a matrix containing the location of each vertex;
2. an matrix describing the element connectivity (the vertices belonging to each element).

For example, suppose `vtx` is a number-of-vertices x 3 matrix of vertices, and `ele` is a number-of-elements x 4 matrix defining the element connectivity, then we may construct a mesh of tetrahedral elements as follows.

```
mesh = Mesh(vtx, ele)
```

The following tabulates some important methods which operate on meshes, see inline help for details of the arguments and return types.

| Function    | Purpose                                                                             |
|-------------|-------------------------------------------------------------------------------------|
| load        | Load a mesh from disk                                                               |
| save        | Save a mesh to disk                                                                 |
| nodecount   | The number of nodes in the mesh                                                     |
| elemcount   | The number of elements in the mesh                                                  |
| dimensions  | The dimensionality of the mesh                                                      |
| boundingbox | Dimensions of a bounding box (square, cube) that encloses the mesh                  |
| fullsize    | The area or volume of the mesh                                                      |
| data        | Returns the vertices, element connectivity and element type of the mesh             |
| surface     | Extract and return the vertices and element connectivity of the surface of the mesh |
| assemble(!) | Build a sparse matrix or a vector of assembled element integrals                    |

### Assembly

TOAST++, like most general purpose finite element libraries, allows the user to assemble system matrices and right hand side vectors by computing local integrals over elements, and assembling these contributions into the global system matrix.

At present, the libTOAST.jl interface only exposes the ability to compute the fully assembled matrices in a single step. One advantage of this approach is that the back-end library can perform multi-threaded assembly (Julia's multi-threading capabilities are still experimental).

The nature of the assembly, and the domain of the integration, is specified by an integration symbol that is defined as one of three enumerations: LinearIntegrals, BilinearIntegrals, and BilinearParamIntegrals. Assembly is performed by one of the following methods:

1. `assemble(mesh, integral [, parameter])`;
2. `assemble!(sysmat, integral [, parameter])`.

Method (1) creates a new `SystemMatrix` object and adds the contributions from the specified integral. Method (2) adds the contribution to an existing `SystemMatrix` object, in-place.

The purpose of the `SystemMatrix` type is described in part (4) of the first example in this document.

#### Assembly types

libTOAST.jl supports the assembly of the following bilinear forms:

| Symbol | Integral                                | Domain |
|--------|-----------------------------------------|--------|
| FF     | Sᵢⱼ = ∫ uᵢ(r) uⱼ(r) dr                  | Ω      |
| DD     | Sᵢⱼ = ∫ ∇uᵢ(r)⋅∇uⱼ(r) dr                | Ω      |
| PFF    | Sᵢⱼ = ∑ᵥ fᵥ ∫ uᵥ(r) uᵢ(r) uⱼ(r) dr      | Ω      |
| PDD    | Sᵢⱼ = ∑ᵥ fᵥ ∫ uᵥ(r) ∇uᵢ(r)⋅∇uⱼ(r) dr    | Ω      |
| BNDPFF | Sᵢⱼ = ∫ uᵢ(r) uⱼ(r) dr                  | δΩ     |
| BNDPFF | Sᵢⱼ = ∑ᵥ fᵥ ∫ uᵥ(r) ∇uᵢ(r)⋅∇uⱼ(r) dr    | δΩ     |

and the following linear forms:

| Symbol | Integral                                | Domain |
|--------|-----------------------------------------|--------|
| F      | Qᵢ = ∫ uᵢ(r) dr                         | Ω      |
| PF     | Qᵢ = ∑ᵥ fᵥ ∫ uᵥ(r) uᵢ(r) dr             | Ω      |
| BNDPF  | Qᵢ = ∑ᵥ fᵥ ∫ uᵥ(r) uᵢ(r) dr             | δΩ     |

### Rasters

TOAST++ is designed to allow solution of the inverse problem in DOT, which consists of estimating the parameters of the PDE from knowledge of solutions. It is convenient to perform this image reconstruction process in a pixel- or voxel-wise representation of the domain, rather than directly in the mesh basis.

To this end, TOAST++ [provides a technique](http://electronicimaging.spiedigitallibrary.org/article.aspx?articleid=1098083) by which to define a raster of pixels/voxels (and indeed, other functions such as radial basis functions) over the support of the underlying domain.

A raster provides three bases, and libTOAST.jl expresses these in three subtypes of an AbstractArray.

| Symbol | Basis and purpose                                                      | Coefficient Type  |
|--------|------------------------------------------------------------------------|-------------------|
| b      | a uniform rasterisation over the bounding box of the mesh              | RasterCoeff       |
| s      | as (b), excludes raster elements which not in the support of the mesh  | SolutionCoeff     |
| g      | as (b), but of higher resolution to improve mapping between bases      | IntermediateCoeff |

To construct a raster of pixels over a given `mesh`

```
julia> raster = Raster(mesh, Pixel)
```

One may then map a function `fm` defined on the mesh (as a `NodalCoeff` type) to the various raster bases

```
julia> fr = RasterCoeff(raster, fm)
julia> fs = SolutionCoeff(raster, fm)
julia> fc = IntermediateCoeff(raster, fm)
```

and back again

```
julia> fm = NodalCoeff(fr)
```

Alternatively, one may map in-place

```
julia> fr = RasterCoeff(raster)
julia> map!(fr, fm)
julia> map!(fm, fr)
```

### Regularisation

TOAST++ provides the functionality to calculate the value, gradient and Hessian of various regularisation functionals. The libTOAST.jl interface provides access to a subset of these functionals. Note that all regularisation methods operate on `SolutionCoeff` types.

| Symbol | Type                    | Value            |
|--------|-------------------------|------------------|
| TK0    | Tikhonov (zeroth-order) | `|| x-x₀ ||²`    |
| TK1    | Tikhonov (first-order)  | `|| ∇(x-x₀) ||²` |
| TV     | Soft Total-Variation    |                  |
| PM     | Perona-Malik            |                  |

The functionals are initialised with a baseline parameter in the solution basis of a raster, e.g.,

```
julia> reg = Regularisation(TK0, x0)
```

subsequently, one may compute the value, gradient and Hessian using the provided methods.

| Function | Purpose                                   |
|----------|-------------------------------------------|
| value    | Value of the regularisation functional    |
| grad     | Gradient of the regularisation functional |
| hess     | Hessian of the regularisation functional  |


# References

[TOAST++ GitHub repository](https://github.com/toastpp/toastpp)

[The Toast++ software suite for forward and inverse modeling in optical tomography](http://biomedicaloptics.spiedigitallibrary.org/article.aspx?articleid=1867204)
