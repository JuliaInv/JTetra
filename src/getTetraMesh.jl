export getTetraMesh

function getTetraMesh(P::Array{Float64},T::Array{Int64})

	T = sort(T,2)

	(edges,faces) = getTetraMeshNumbering(T)

	nc = size(T,1)
  nf = size(faces,1)
  ne = size(edges,1)
	nn = size(P,1)

	Div   = spzeros(0,0)
	Curl  = spzeros(0,0)
	Grad  = spzeros(0,0)
  
  V = Float64[]
  F = Float64[]
  L = Float64[]

	return TetraMesh(P,T,nc,nf,ne,nn,Div,Curl,Grad,V,F,L,faces,edges)
  
end

function getTetraMesh(PFile::String,TFile::String)

  file  = open(PFile, "r")
  lines = readlines(file)
  close(file)

  n = length(lines)

  P = zeros(Float64, n, 3)

  for i = 1:n
    line = split(lines[i])
    P[i,1] = parse(Float64, line[1])
    P[i,2] = parse(Float64, line[2])
    P[i,3] = parse(Float64, line[3])
  end

  file  = open(TFile, "r")
  lines = readlines(file)
  close(file)

  n = length(lines)

  T = zeros(Int64, n, 4)

  for i = 1:n
    line = split(lines[i])
    T[i,1] = parse(Int64, line[1])
    T[i,2] = parse(Int64, line[2])
    T[i,3] = parse(Int64, line[3])
    T[i,4] = parse(Int64, line[4])
  end
  
  return getTetraMesh(P,T)

end