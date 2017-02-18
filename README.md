# TOAST.jl

TOAST.jl is a low-level inteface to the [TOAST++](https://github.com/toastpp/toastpp) library. Whilst TOAST++ is designed as an end-to-end solution for forward modeling and image reconstruction in Diffuse Optical Tomography (DOT), the TOAST.jl interface simply provides an interface to its underlying  finite-elemenent, raster mapping, and regularisation functionality. DOT specific functionality is implemented by a seperate Julia package.

## Compatibilty

TOAST.jl directly calls the underlying C++ library using [Cxx.jl](https://github.com/Keno/Cxx.jl) which is currently supported only on Linux and OS X. Windows supported is expected with the release of Julia 0.6.

## Installation

```
# Dependencies
Pkg.add("Cxx")

# TOAST.jl
Pkg.clone("git@github.com:samuelpowell/TOAST.jl.git")
Pkg.build("TOAST")
```

The TOAST.jl build process will download the requisite binaries and the header files required to call the library, these are stored in the package directory. You can test that TOAST.jl is properly configured by running the unit tests.

```
# Run all unit-tests
Pkg.test("TOAST")
```

## A first example

By means of an introduction, we will use TOAST.jl to solve a lossy diffusion equation with Neumann boundary conditions.

### 1. Load the TOAST module

To begin, load the TOAST module.

```
julia> using TOAST
Default thread count set to 8
```

The TOAST++ backend reports that it has initialised a pool of eight threads in which it will perform tasks such as system matrix assembly.

### 2. Define the domain

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

There are two parameters in the problem, μa, and μs. In TOAST++ it is assumed that functions (be they parameters, or source distributions) can be provided as a set of nodal coefficients in the same basis as the shape functions used for the finite-element solution. That is to say that a function `f(r)` defined over the domain is represented as

```
f(r) = ∑_k f_k u_k(r)
```

Functions in nodal basis are represented in `NodalCoeff` types, which hold a reference to the `mesh`.

```
julia> μ = fill(NodalCoeff, mesh, 0.01);
julia> κ = fill(NodalCoeff, mesh, 0.3);
julia> ζ = fill(NodalCoeff, mesh, 0.5);
```

Evidently the length of the nodal coefficient vectors is equal to the number of nodes in the mesh.

```
julia> nodecount(mesh), length(μa)
(3511,3511)
```

### 4. Build the system matrix

For the specified problem the system matrix is given by

```
S = K + M + F
```
where

```
Kᵢⱼ = ∫_Ω κ(r) ∇uᵢ(r)⋅∇uⱼ(r) dr = ∑_k κ_k ∫_Ω u_k(r) ∇uᵢ(r)⋅∇uⱼ(r) dr
Mᵢⱼ = ∫_Ω μ(r) uᵢ(r) uⱼ(r) dr = ∑_k μ_k ∫_Ω u_k(r) uᵢ(r) uⱼ(r) dr
Fᵢⱼ = 1/(2A) ∫_δΩ uᵢ(r) uⱼ(r) dr
```

To build these matrices

```
julia> K = assemble(mesh, PDD, κ);
julia> M = assemble(mesh, PFF, μ);
julia> F = assemble(mesh, BNDPFF, ζ);
```

The first input to `assemble` is the `mesh` over which the integrals are performed, the second defined the type of integral, and the third is the parameter. The available integrals are detailed later, but for now it will suffice to understand that `P` represents a paraemter in the mesh basis, `F` represents a basis function, `D` represents a derivative of a basis function, and the `BND` prefix indicates the domain of the integration is only over the boundary.

The result of each `assembly` operation is a `TOAST.SystemMatrix` type, which is a wrapper around a compressed sparse row (CSR) matrix used by TOAST. Since Julia's native sparse arrays are compressed sparse column (CSC) TOAST.jl uses a thin wrapper around the native format until the final matrix is ready for use. This allows us to instead do the following.

```
julia> S = SystemMatrix(mesh)
julia> assemble!(S, PDD, κ);
julia> assemble!(S, PFF, μ);
julia> assemble!(S, BNDPFF, ζ);
```

Here only a single sparse matrix is allocated, and the results of the assembly operations are added in-place. Note that these mutating methods do not require the mesh to be specified, as a reference is held in the `SystemMatrix` type. When we are finally ready to use the system matrix in Julia, we can do so as follows.

```
sysmat = sparse(S)
```

### 4. Build a source term

Equipped with a system matrix encoding the underlying PDE in a finite-dimensional subspace, we can solve the system against a source term (right hand side). The general process is the same as assemblying a system matrix, however we specify linear integrals. 

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
julia> q = assemble(mesh, PF, source)
```

### 5. Solve the system

The resultant system `Sϕ = q` can be solved using your method of preference. For example, we can simply rely upon Julia's backslash operator to perform a sensible direct solve

```
julia> ϕ = S\q
```

et voilà.

## Concepts

## Types and Methods

### Meshes and the FE subspace

In TOAST the mesh and the underlying finite-element subspace are tightly coupled. 

### Assembly

#### Assembly types

### Rasters

### Regularisation

