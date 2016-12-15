module JTetra

importall jInv.Mesh
using jInv.Utils

export TetraMesh

type TetraMesh <: AbstractMesh

	Points::Array{Float64}      # Points of tesselation
  Tetras::Array{Int64}        # Connectivity map

	nc::Int  # number of cells
	nf::Int  # number of faces
	ne::Int  # number of edges
	nn::Int  # number of nodes

	Div::SparseMatrixCSC
	Curl::SparseMatrixCSC
	Grad::SparseMatrixCSC

	V::Array{Float64,1} # cell volume
	F::Array{Float64,1} # face area
	L::Array{Float64,1} # edge lengths

	faces::Array # FaceNumbering
	edges::Array # EdgeNumbering

end # type TetraMesh

include("getTetraMesh.jl")
include("getTetraMeshNumbering.jl")
include("findRow.jl")
include("getNodalMassMatrix.jl")
include("getEdgeMassMatrix.jl")
include("getFaceMassMatrix.jl")
include("getNodalGradientMatrix.jl")
include("getCurlMatrix.jl")
include("getDivergenceMatrix.jl")
include("getSizes.jl")
include("getEdgeIntegralOfPolygonalChain.jl")

import Base.clear!
function clear!(M::TetraMesh)
	M.Div  = clear!(M.Div)
	M.Curl = clear!(M.Curl)
	M.Grad = clear!(M.Grad)
end # function clear

export clear!

end
