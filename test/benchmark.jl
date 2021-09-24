
result = [6.83 205 6.83 6 4; 204 4 4.83 9.83 9; 2.83 2 7 11.9 17; 2 0 4 207 221; 2.83 NaN 5.66 209 413]
points = [0.0 0 0 0 2; 0 0 0 0 0; 0 0 0 0 0; 0 1 0 0 0; 0 0 0 0 0]
initial = [8.0 8 8 8 4; 8 8 8 8 8; 8 8 8 8 8; 0 0 8 8 8; 0 0 8 8 8]
friction = [1.0 200 1 1 1; 200 1 1 4 4; 1 1 4 4 4; 1 1 3 200 200; 1 Inf 3 200 4]

lpoints = repeat(points, 10, 10)
linitial = repeat(initial, 10, 10)
lfriction = repeat(friction, 10, 10)

llpoints = repeat(points, 100, 100)
llinitial = repeat(initial, 100, 100)
llfriction = repeat(friction, 100, 100)

lllpoints = repeat(points, 1000, 1000)
lllinitial = repeat(initial, 1000, 1000)
lllfriction = repeat(friction, 1000, 1000)

# @test result == spread(points, initial, friction)
spread(points, initial, friction, res=2);


x = [5, 2, 6, 8] .* 100
y = [2, 6, 7, 8] .* 100
points = zeros(1000, 1000)
points[CartesianIndex.(x, y)] .= 1

friction = ones(1000, 1000)
initial = zeros(1000, 1000)
friction[580:600, begin:450] .= Inf

spread2(points, initial, friction)
