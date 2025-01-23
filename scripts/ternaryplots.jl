using CairoMakie, AnophelesSDMs

convert_to_2d(r1)

fig = Figure()
ax = Axis(fig[1,1])

text!(ax, convert_to_2d([0.4, 0.6, 0]); text = "     " * string(x), align = (:center, :center), rotation = pi/3 * rot, color = :red)

ticks = string.(linerange)

for pos in 1:3
    tick_positions = [[l, -0.05, 1 - l] for l in linerange]
    text!(ax, convert_to_2d.(tick_positions); text = ticks, align = (:left, :center))
end

fig

lines!(ax, lines, color = :black)


points = [[0.1, 0.1, 0.8], [0.5, 0.4, 0.1]]
scatter!(ax, Tuple.(convert_to_2d.(points)))
fig
[R * point for point in points]