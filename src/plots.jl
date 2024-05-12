using Plots

# Figure 4 Panel A

minp = [minimum(p[i, Sx[i, :] .> 0]) for i in 1:N]

function fig4a()
    plot(y, minp, color=:red, linewidth=2, label="viability threshold")
    plot!(y, y, linewidth=2, label="aggregate shock")

    for j in 1:10:M
        if j == 1
            # Add a label only on the first iteration
            plot!(y, p[:, j], linestyle=:dash, linewidth=1.2, label="productivity by type")
        else
            # No label for other lines
            plot!(y, p[:, j], linestyle=:dash, linewidth=1.2, label="")
        end
    end


    xlabel!("Aggregate shock")
    ylabel!("Match productivity")
    plot!(xlims=(minimum(y), maximum(y)), ylims=(minimum(0.65), maximum(1.8)))
 
    plt = plot!()

    savefig(plt, "figures/fig4a_plot.png")

end

fig4a()


# Figure 4 Panel b

function fig4b()
    
    wFb = ones(N, 1)

    ability_values = B * x .+ C
    
    unemployed_proportion = ((wFb' * ux)' .* l) ./ (wFb' * ux * l)

    # Initialize plot with unemployed proportion
    plot(ability_values, unemployed_proportion, linestyle=:dashdot, color=:blue, linewidth=2, label="Unemployed")

    plot!(ability_values, l, color=:red, linewidth=2, label="All")


    xlabel!("Ability")
    ylabel!("Proportion")
    
    plot!(xlims=(minimum(0.7), maximum(1.7)), ylims=(minimum(0), maximum(0.02)))

    # Save the figure
    plt = plot!()  
    savefig(plt, "figures/fig4b_plot.png")  
end

# Call the function to display and save the plot
fig4b()


# figure 6

function fig6()
    
    # Plotting yt vs. ut with stars
    plot(yt, ut, seriestype=:scatter, shape=:star, color=:blue, legend = false)

    # Adding y vs. u with a red line
    plot!(y, u, color=:red, linewidth=2)

    # Adding labels
    ylabel!("Unemployment Rate")
    xlabel!("Aggregate Shock")

    # Display the plot
    plt = plot!()  # Make sure to collect the plot object

    # Save the plot
    savefig(plt, "figures/fig6_plot.png")  # Saves as EPS file
end

# Call the function to execute the plotting and saving
fig6()


# Figure 10
function fig10()
    
    # Calculate the cumulative sum of 'l'
    Fl = cumsum(l)

    # Define the quantiles
    q = [0.1, 0.25, 0.5, 0.75, 0.9]

    # Initialize 
    mqua = zeros(Int, length(q))


    # Determine the indices corresponding to each quantile
    for i in 1:length(q)
        mqua[i] = findlast(Fl .<= q[i] * last(Fl))
    end

    # get min and max wages 
    wmin = wd[:wmin]
    wmax = wd[:wmax]


    # Initialize the plot
    plt = plot(legend = false)

    # Loop through different groups
    for m in 1:5
        # Plot minimum wages with circle markers
        plot!(y, wmin[:, mqua[m]], color=:red, line=:solid, marker=:circle, 
              markersize=2 * m, label="Min Wage Group $m")

        # Plot maximum wages with square markers
        plot!(y, wmax[:, mqua[m]], color=:blue, line=:solid, marker=:square, 
              markersize=2 * m, label="Max Wage Group $m")
    end

    # Adding labels and customize axes
    xlabel!("Productivity")
    ylabel!("Wages")

    # Tighten the axis around the data
    plot!(xlims=(minimum(0.93), maximum(1.08)), ylims=(minimum(0.5), maximum(1.4)))

    # Display the plot
    plt = plot!()

    # Save the plot
    savefig(plt, "figures/fig10_plot.png") 

end

# Call the function to generate and save the plot
fig10()