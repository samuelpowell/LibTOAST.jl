# TOAST.jl

TOAST.jl is a low-level inteface to the [TOAST++](https://github.com/toastpp/toastpp) library. Whilst TOAST++ is designed as an end-to-end solution for forward modeling and image reconstruction in Diffuse Optical Tomography (DOT), the TOAST.jl interface simply provides an interface to its underlying  finite-elemenent, raster mapping, and regularisation functionality. DOT specific functionality is implemented in a seperate Julia package.

## Compatibilty

TOAST.jl directly calls the underlying C++ library using [Cxx.jl](https://github.com/Keno/Cxx.jl) which is currently supported only on Linux and OS X. Windows support is expected with the release of Julia 0.6.

## Installation

Ensure that you have a working Cxx.jl installation, then install and build TOAST.jl package as follows.

```
julia> Pkg.clone("git@github.com:samuelpowell/TOAST.jl.git")
julia> Pkg.build("TOAST")
```

The TOAST.jl build process will download the requisite binaries and the header files required to call the library, these are stored in the package directory. You can test that TOAST.jl is properly configured by running the unit tests.

```
Pkg.test("TOAST")
```

## A first example

By means of an introduction, we will use TOAST.jl to solve a lossy diffusion equation over a domain Ω

```
[-∇⋅κ(r)∇ + μ(r)] ϕ(r) = 0 (r ∈ Ω),
```

with a Robin condition on the boundary

```
ϕ(r) + γ n⋅∇ϕ(r) = q(r) (r ∈ δΩ),
```

where `n` is the unit outward normal to the boundary. To solve this PDE numerically by the finite element method we subdivide the domain into a set of non-overlapping elements joined at N vertex nodes, and represent the discrete solution in a piecewise linear basis

```
ϕ(r) ≈ ϕʰ(r) = ∑_i ϕᵢ uᵢ(r)
```

We represent approximations to the parameters γ(r), κ(r), and μ(r) in the same basis. Application of the Galerkin method to the continuous PDE yields the following matrix equation,

```
S ϕ = Q
```

where the system matrix `S = K + M + F` is calculated as

```
Kᵢⱼ = ∑_k κ_k ∫_Ω  u_k(r) ∇uᵢ(r)⋅∇uⱼ(r) dr
Mᵢⱼ = ∑_k μ_k ∫_Ω  u_k(r)  uᵢ(r) uⱼ(r) dr
Fᵢⱼ = ∑_k γ_k ∫_δΩ u_k(r) uᵢ(r) uⱼ(r) dr
```
and

```
Qᵢ = ∑_i qᵢ ∫_δΩ uᵢ(r)
```

We will now proceed to find the vector of nodal solution coefficeints, ϕ, by assembling the relevant matrices and source vectors using the TOAST.jl interface to TOAST++.

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
julia> meshdir = joinpath(TOAST._jl_toast_dir, "test", "2D", "meshes")

julia> mesh = Mesh(joinpath(meshdir, "circle25_32.msh"));
INFO: Loaded 250114 bytes from /Users/spowell/.julia/v0.5/TOAST/deps/toastpp-2.0.2/test/2D/meshes/circle25_32.msh.
INFO: Initiliasing mesh.
INFO: Calculating sparsity pattern.
INFO: Number of vertices: 3511, number of nonzeros: 24211
```

TOAST.jl reports that it has loaded the mesh. In TOAST++ the mesh and the finite-element subspace are tightly coupled: a description of the mesh includes the form of the finite-element shape functions; as such TOAST.jl is able to precompute the sparsity pattern which will be used for system matrix assembly.

### 3. Define the parameters

There are three parameters in the problem, κ(r), μ(r), and ζ(r). We directly provide the discrete form of the parameters as a set of nodal coefficients in the same basis as the shape functions used for the finite-element solution. 

Functions in nodal basis are represented in `NodalCoeff` types, which hold a reference to the `mesh`.

```
julia> μ = fill(NodalCoeff, mesh, 0.01);
julia> κ = fill(NodalCoeff, mesh, 0.3);
julia> ζ = fill(NodalCoeff, mesh, 0.5);
```

### 4. Build the system matrix

To build the system matrix we first intatiate a `SystemMatrix`

```
julia> S = SystemMatrix(mesh)
```

which is a thin wrapper around a compressed sparse row (CSR) matrix used by TOAST. Since Julia's native sparse arrays are compressed sparse column (CSC) TOAST.jl uses this wrapper until the final matrix is ready for use, this allows us to add components of the system matrix in-place without additional allocations.

We now assemble and add each of the components to the system matrix.

```
julia> assemble!(S, PDD, κ);
julia> assemble!(S, PFF, μ);
julia> assemble!(S, BNDPFF, ζ);
```

The second input to the `assemble!` method defines the integral to be performed. The available integrals are enumerated later, but for now it will suffice to understand that `P` represents a paraemter in the mesh basis, `F` represents a basis function, `D` represents a derivative of a basis function, and the `BND` prefix indicates the integration is only over the boundary.

When we are finally ready to use the system matrix in Julia, we can do so as follows.

```
sysmat = sparse(S)
```

### 4. Build a source term

Equipped with a system matrix describing the discrete form of the PDE we can solve the system against an arbitrary source term (right hand side). Building the source term follows a similar process to assembling the system matrix, though there are some additional steps since the source is spatially varying.

First, we construct the source function in space. To begin, we need to extract the position of the vertices in the mesh

```
julia> vtx, = data(mesh)
```

where vtx will be a `nodecount(mesh) x dimension(mesh)` matrix of vertices. We may now construct a source function which depends upon this geometry

```
julia> source = NodalCoeff(mesh, 1./((vtx[:,1]^2 + vtx[:,2]^2) + 0.1 )
```

before projecting this into the mesh

```
julia> Q = assemble(mesh, P, source)
```

### 5. Solve the system

The resultant system `Sϕ = Q` can be solved using your method of preference. For example, we can simply rely upon Julia's backslash operator to determine the properties of the system matrix and undertake a suitable direct solve.

```
julia> ϕ = sysmat\q
```

## Concepts

## Types and Methods

### Meshes and the FE subspace

In TOAST the mesh and the underlying finite-element subspace are tightly coupled. 

### Assembly

#### Assembly types

TOAST.jl supports the assembly of the following bilinear forms

* FF: Sᵢⱼ = ∫_Ω  uᵢ(r) uⱼ(r) dr
* DD: Sᵢⱼ = ∫_Ω  ∇uᵢ(r)⋅∇uⱼ(r) dr
* PFF: Sᵢⱼ =  ∑_k f_k ∫_Ω  u_k(r) uᵢ(r) uⱼ(r) dr
* PDD:  Sᵢⱼ = ∑_k f_k ∫_Ω  u_k(r) ∇uᵢ(r)⋅∇uⱼ(r) dr
* BNDFF: Sᵢⱼ = ∫_δΩ   uᵢ(r) uⱼ(r) dr
* BNFPFF: Sᵢⱼ = ∑_k f_k ∫_δΩ  u_k(r) ∇uᵢ(r)⋅∇uⱼ(r) dr

and the following linear forms

* F: Qᵢ = ∫_Ω  uᵢ(r) dr
* PF: Qᵢ =  ∑_k f_k ∫_Ω  u_k(r) uᵢ(r) dr
* BNDF: Qᵢ = ∑_k f_k ∫_δΩ  u_k(r) uᵢ(r) dr

### Rasters

TOAST++ is designed to allow solution of the inverse problem in DOT, which consists of estimating the parameters of the PDE from knowledge of solutions. It is convenient to perform this reconstruction in a pixel- or voxel-wise representation of the domain, rahter than directly in the mesh basis.

To this end, TOAST++ provides a technique by which to define a raster of pixels/voxels (and indeed, other functions such as radial basis functions) over support of the underlying domain.

A raster provides three bases:

1. raster basis (b), which is a uniform rasterisation over the bounding box of the mesh
2. solution basis (s), which is the same as the raster, but excludes raster elements which are not in the support of the mesh
3. intermediate basis (g), which is a higher resolution version of the raster basis, used to improve the quality of the mapping between.


To construct a raster of pixsels over a given `mesh`

```
raster = Raster(mesh, Pixel)
```

One may then map a function `fm` defined on the mesh (as a NodalCoeff type) to the various raster bases

```
fr = RasterCoeff(raster, fm)
fs = SolutionCoeff(raster, function)
fc = IntermediateCoeff(raster, function)
```

and back again

```
fm = NodalCoeff(fr)
```

Alternatively, one may map in-place

```
fr = RasterCoeff(raster)
map!(fr, fm)
map!(fm, fr)
```

### Regularisation

TOAST++ provides the functionality to calculate the value, gradient and Hessian of various regularisation functionals. The TOAST.jl interface provides access to a subset of these functionals:

* TK0: zeroth-order Tikhonov
* TK1: first-order Tikhonov
* TV: total-variation (soft)
* PM: Perona-Malik

The functionals are initialised with a baseline parameter in the solution basis of a raster, e.g.,

```
reg = Regularisation(TK0, x0)
```

subsequently, one may compute the value, gradient and Hessian as follows

```
value(reg, x)
grad(reg, x)
hess(reg, x)
```