using Plots

# Define the x range for two periods of cosine (0 to 4π)
x = 0:0.01:6π
y = cos.(x)

# Plot the cosine curve
plot(x, y, label="cos(x)", linewidth=2, legend=:topright)

# Add 8 equally spaced points over the domain
n_marks = 15
x_marks = range(0, stop=6π, length=n_marks)
y_marks = cos.(x_marks)

plt1 = Plots.plot(x,y, color = :blue, label = "", lw = 3, framestyle = :box, yticks = false,tickfontsize = 14,guidefontsize = 18)
xticks!(plt1, (0:2π:6π,["0","a","2a","3a"]))
xlabel!(plt1, "x")
# Add small horizontal marks (line segments) on top of the cosine
for (xm, ym) in zip(x_marks, y_marks)
    plot!(plt1, [xm - 0.2, xm + 0.2], [ym, ym], color=:red, linewidth=4, label = "")
end

# Final plot display
display(plt1)