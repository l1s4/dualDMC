# Input documentation:
# named_values: a named numeric vector
# sigma_old, sigma_new: the previous and target diffusion constants
# t_from_to: scaling of time (options: ms->s, s->ms, or none)
convert_prms <- function(
  named_values,
  sigma_old = 4,
  sigma_new = 1,
  t_from_to = "ms->s"
) {
  # Some rough input checks
  stopifnot(is.numeric(named_values), is.character(names(named_values)))
  stopifnot(is.numeric(sigma_old), is.numeric(sigma_new))
  t_from_to <- match.arg(t_from_to, choices = c("ms->s", "s->ms", "none"))

  # Internal conversion function (takes a name and value pair, and transforms it)
  internal <- function(name, value) {
    name <- match.arg(
      name,
      choices = c("muc", "b", "non_dec", "sd_non_dec", "tau", "a", "A", "alpha")
    )

    # 1. scale for the diffusion constant
    if (name %in% c("muc", "b", "A")) {
      value <- value * (sigma_new / sigma_old)
    }

    # 2. scale for the time
    # determine the scaling per parameter (assuming s->ms)
    scale <- 1
    if (name %in% c("non_dec", "sd_non_dec", "tau")) {
      scale <- 1000
    }
    if (name %in% c("b", "A")) {
      scale <- sqrt(1000)
    }
    if (name %in% c("muc")) {
      scale <- sqrt(1000) / 1000
    }

    # adapt, depending on the t_from_to argument
    if (t_from_to == "ms->s") {
      scale <- 1 / scale
    }
    if (t_from_to == "none") {
      scale <- 1
    }

    value <- value * scale
  }

  # Apply the internal function to each element
  converted_values <- mapply(FUN = internal, names(named_values), named_values)

  return(converted_values)
}


