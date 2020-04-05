library(ggvis)


head(mtcars)

#basic scatterplot
ggvis(data = mtcars, ~mpg, ~wt, fill := "red")
ggvis(data = mtcars, ~mpg, ~wt, fill = ~cyl)

layer_points(ggvis(mtcars, ~mpg, ~wt))

mtcars %>% ggvis(~mpg, ~wt) %>% layer_points()

mtcars %>%
  ggvis(~mpg, ~wt) %>%
  layer_points() %>%
  layer_smooths()

mtcars %>%
  ggvis(~mpg, ~wt) %>%
  layer_points() %>%
  layer_smooths(se=TRUE)

mtcars %>%
  ggvis(~mpg, ~wt) %>%
  layer_points(fill= ~factor(cyl))

mtcars %>%
  ggvis(~mpg, ~wt) %>%
  layer_points(fill= ~factor(cyl)) %>% 
  group_by(cyl) %>%
  layer_model_predictions(model = "lm", se=TRUE)

mtcars %>%
  ggvis(~mpg, ~wt) %>%
  layer_smooths(span = input_slider(0.5, 1, value = 1)) %>%
  layer_points(size := input_slider(100,1000, value=100), fill= ~factor(cyl))


mtcars %>%
  ggvis(x = ~wt) %>%
  layer_densities(
    adjust = input_slider(.1, 2, value = 1, step = .1, label = "Bandwidth adjustment"),
    kernel = input_select(
      c("Gaussian" = "gaussian",
        "Epanechnikov" = "epanechnikov",
        "Rectangular" = "rectangular",
        "Triangular" = "triangular",
        "Biweight" = "biweight",
        "Cosine" = "cosine",
        "Optcosine" = "optcosine"),
      label = "Kernel")
  )