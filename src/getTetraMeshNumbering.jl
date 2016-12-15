function getTetraMeshNumbering(Tetras)

# tetrahedron connectivity
t = Tetras

nc = size(t,1)
# global edge vertex numbers
edges = zeros(Int64,6*nc,2)

cnt = 1
for i=1:nc
	edges[cnt,1] = t[i,1]
	edges[cnt,2] = t[i,2]; cnt += 1

	edges[cnt,1] = t[i,1]
	edges[cnt,2] = t[i,3]; cnt += 1

	edges[cnt,1] = t[i,1]
	edges[cnt,2] = t[i,4]; cnt += 1

	edges[cnt,1] = t[i,2]
	edges[cnt,2] = t[i,3]; cnt += 1

	edges[cnt,1] = t[i,2]
	edges[cnt,2] = t[i,4]; cnt += 1

	edges[cnt,1] = t[i,3]
	edges[cnt,2] = t[i,4]; cnt += 1
end


# collate duplicate edges
edges = unique(sortrows(edges), 1)

# face degrees of freedom
faces = zeros(Int64,4*nc,3)
cnt = 1
for i=1:nc
	faces[cnt,1] = t[i,1]
	faces[cnt,2] = t[i,2]
	faces[cnt,3] = t[i,3]; cnt += 1

	faces[cnt,1] = t[i,1]
	faces[cnt,2] = t[i,2]
	faces[cnt,3] = t[i,4]; cnt += 1

	faces[cnt,1] = t[i,1]
	faces[cnt,2] = t[i,3]
	faces[cnt,3] = t[i,4]; cnt += 1

	faces[cnt,1] = t[i,2]
	faces[cnt,2] = t[i,3]
	faces[cnt,3] = t[i,4]; cnt += 1
end

# global face vertex numbers
faces = unique(sortrows(faces), 1)

return edges, faces

end
