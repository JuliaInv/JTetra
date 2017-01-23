using JTetra
using Base.Test

println("==== Test persistency of linear operators ====")

# read tetrahedral mesh
mesh = getTetraMesh("points.txt","tetras.txt")

# scalar and tensor cell properties for mass matrices
a = ones(mesh.nc)
A = [ones(mesh.nc, 3) zeros(mesh.nc, 3)]

# build all matrices
Grad = getNodalGradientMatrix(mesh)
Curl = getCurlMatrix(mesh)
Div  = getDivergenceMatrix(mesh)

L = getEdgeLength(mesh)
F = getFaceArea(mesh)
V = getCellVolume(mesh)

Mn  = getNodalMassMatrix(mesh, a)
Mea = getEdgeMassMatrix(mesh, a)
MeA = getEdgeMassMatrix(mesh, A)
Mfa = getFaceMassMatrix(mesh, a)
MfA = getFaceMassMatrix(mesh, A)

@test Mea == MeA
@test Mfa == MfA

# now all fields should be non-empty
@test !isempty(mesh.Grad)
@test !isempty(mesh.Curl)
@test !isempty(mesh.Div)
@test !isempty(mesh.L)
@test !isempty(mesh.F)
@test !isempty(mesh.V)

clear!(mesh)

# now all three matrices should be empty
@test isempty(mesh.Grad)
@test isempty(mesh.Curl)
@test isempty(mesh.Div)

print("Passed!\n")


println("=== Test getEdgeIntegralOfPolygonalChain ===")

# define one transmitter
trx = [
  -5.0 0.0 0.0
   5.0 0.0 0.0]

Src = getEdgeIntegralOfPolygonalChain(mesh, trx, normalize = false)

# define a receiver coil - closed polygon
rcv1 = [
  -5.0 -5.0 0.0
  -5.0  5.0 0.0
   5.0  5.0 0.0
   5.0 -5.0 0.0
  -5.0 -5.0 0.0]

Obs1 = getEdgeIntegralOfPolygonalChain(mesh, rcv1, normalize = true)

# define another receiver - open polygon
rcv2 = [
  -5.0 -5.0 0.0
  -5.0  5.0 0.0]

Obs2 = getEdgeIntegralOfPolygonalChain(mesh, rcv2, normalize = true)

print("Passed!\n")
