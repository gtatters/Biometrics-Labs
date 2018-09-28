# Base graphics
# Base graphics are very flexible and allow a great deal of customisation, 
# with many individual functions available. However, they lack a coherent 
# underlying framework and, for visualizing highly structured data, are 
# outclassed by lattice and ggplot2.

#Quick reference info:
demo("graphics") # Demonstration of graphics in R
?plot            # Help page for main plot function

# Useful plotting functions: 
#  lines, points, abline, curve, text, rug, legend
#  segments, arrows, polygon


# Create some data for plotting: 
x <- 10 + (1:20)/10
y <- x^2 + rnorm(length(x))     # Add Gaussian random number
plot(x, y)
curve(x^2, add=TRUE, lty=2)     # Add dashed line showing y=x^2
plot(x, y, type="l", col="blue")  # Plot as blue line (try 'type="o"')
plot(x, y, type="l", log="xy")   # Plot as line with log X & Y axes:
abline(v=11, lty=3)    # Add vertical dotted line
text(11.5, 120, "Hello")   # Add annotation
legend("topleft", inset=0.05, "data", pch=1, col="blue", bty="n")  # Add a legend
#Different point styles: 
  plot(x, y, pch=2, col="red")  # Hollow triangles
plot(1:10, rep(1, 10), pch=LETTERS)  # Can also use any character
example("pch")  # Show point style examples

# Plot symbols and colours can be specified as vectors, to allow individual specification for each point. R uses recycling of vectors in this situation to determine the attributes for each point, i.e. if the length of the vector is less than the number of points, the vector is repeated and concatenated to match the number required.

# Single plot symbol (see "?points" for more) and colour (type "colours()" or "colors()" for the full list of predefined colours): 
plot(x, y, pch=2, col="red")  # Hollow triangles
plot(x, y, pch=c(3, 20), col=c("red", "blue"))  # Blue dots; red "+" signs
plot(x, y, pch=1:20)   # Different symbol for each point
#Create vector of contiguous colours in a rainbow palette:
  col <- rainbow(length(x))
plot(x, y, col=col)
#Label axes:
  plot(x, y, xlab="Some data", ylab="Wibble")
# Axis limits are controlled by xlim and ylim, which are vectors of the minimum and maximum values, respectively.

# Specify axis limits: 
  plot(x, y, xlim=c(11, 12), ylim=c(0, 150))
