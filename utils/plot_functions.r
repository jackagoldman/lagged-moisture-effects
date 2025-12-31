
mask_too_far <- function(model, x, y, n = 100, dist = 0.05) {
  m <- model.frame(model)

  # check presence
  miss <- setdiff(c(x, y), names(m))
  if (length(miss) > 0) {
    stop("Variable(s) not found in model frame: ", paste(miss, collapse = ", "))
  }

  # helper to coerce to numeric vector, handling matrices/list-columns
  as_num_vec <- function(z) {
    if (is.matrix(z) || is.data.frame(z)) return(as.numeric(z))
    if (is.list(z)) return(unlist(lapply(z, function(el) if (length(el) == 1) as.numeric(el) else as.numeric(el))))
    as.numeric(z)
  }

  x_vals <- as_num_vec(m[[x]])
  y_vals <- as_num_vec(m[[y]])

  # keep only finite observed values
  ok <- is.finite(x_vals) & is.finite(y_vals)
  if (!any(ok)) stop("No finite observed pairs found for variables '", x, "' and '", y, "'. Inspect your data.")
  x_obs <- x_vals[ok]
  y_obs <- y_vals[ok]

  # sequences for grid (now guaranteed finite)
  x_seq <- seq(min(x_obs), max(x_obs), length.out = n)
  y_seq <- seq(min(y_obs), max(y_obs), length.out = n)

  grd <- expand.grid(x = x_seq, y = y_seq)
  names(grd) <- c(x, y)

  # compute mask: TRUE where grid point is too far from any observed pair
  tf_mask <- mgcv::exclude.too.far(grd[[x]], grd[[y]], x_obs, y_obs, dist = dist)

  # return logical vector aligned to grid rows (length = n*n)
  tf_mask
}
#' get lagged smooth estimates from a GAM model
#' @param model, a fitted GAM model
#' @param lag_var, name of the lag variable in the model
#' @param smooth_var, name of the smooth variable in the model
#' @param weights, name of the weights variable in the model
#' @param n, number of points to evaluate the smooth over
#' @return a list with two dataframes: smooth_df (all estimates) and mask (significant estimates)

get_lag_estimates <- function(model, lag_var = "L", smooth_var = "vpd", weights = "weights_btl", n = 100) {
  
  # Create the smooth specification
  smooth = paste0("te(", smooth_var, ",", lag_var,"):",  weights)


  tf <- mask_too_far(model = model, x = lag_var, y = smooth_var, n = 100, dist = 0.1)


  smooth_df <- gratia::smooth_estimates(
    object = model,
    smooth = smooth,
    partial_match = TRUE,
    unconditional = TRUE,
    n = 100
  ) |>
    gratia::add_confint() |> 
# Flag significant locations where CI excludes 0 (null effect)
    tibble::add_column(tf) |> 
    dplyr::mutate(across(c(.estimate, .lower_ci, .upper_ci), ~ if_else(tf, NA_real_, .))) |>
    dplyr::mutate(sig = ifelse(.lower_ci > 0 | .upper_ci < 0, TRUE, FALSE)) 
    
return(smooth_df)
}


#' plot heatmap of lagged smooth effects from a GAM model
#' @param smooth_df, dataframe of smooth estimates with confidence intervals
#' @param mask, dataframe of significant smooth estimates (CI excludes 0)



plot_lag_smooth <- function(model, 
  smooth_df, 
  smooth_var, 
  lag_var, 
  weights) {
  ylabel <- switch(smooth_var,
                   vpd = "Vapour Pressure Deficit (z-score)",
                   cmi = "Climate Moisture Index (z-score)",
                   smooth_var)
  response <- as.character(formula(model)[[2]])
  response_lab <- if (response == "median") "P(median)" else if (response == "pct_90") "P(extreme)" else response
  #if response = "median" set breaks  = 1:30, lmits = c(1,30) but if response = "pct_90" set breaks = 1:5, limits = c(1,5)
  breaks <- if (smooth_var == "vpd") 1:30 else if (smooth_var == "cmi_wy") 1:5 else waiver()
  limits <- if (smooth_var == "vpd") c(1,30) else if (smooth_var == "cmi_wy") c(1,5) else NULL
  # base heatmap mapping columns by name
  p <- ggplot(smooth_df, aes(x = .data[[lag_var]], y = .data[[smooth_var]], fill = .estimate)) +
    geom_tile() +
      scale_fill_viridis_c(
        name = response_lab,
        option = "A",
        guide = guide_colorbar(
          direction = "horizontal",
          title.position = "top",
          title.hjust = 0.5
        )
      ) +
      theme(legend.position = "bottom") +
      scale_x_continuous(
        "Lag (days before fire start)",
        breaks = breaks,
        expand = c(0, 0)
      ) +
    scale_y_continuous(name = ylabel, expand = c(0, 0)) +
    theme_bw() +
    theme(
      panel.spacing = grid::unit(0, "pt"),
      plot.margin = grid::unit(c(0, 0, 0, 0), "pt")
    )

  # outline significant cells (sig == TRUE)
  mask <- dplyr::filter(smooth_df, sig == TRUE)

  if (nrow(mask) > 0) {
     p <- p +
      geom_tile(
        data = mask,
        aes(x = .data[[lag_var]], y = .data[[smooth_var]]),
        color = "black",
        size = 0.75,
        linejoin = "round",
        fill = NA,
        inherit.aes = FALSE
      ) +
      geom_tile(data = mask, aes(x = .data[[lag_var]], y = .data[[smooth_var]], fill = .estimate), inherit.aes = FALSE)
  }
  
  # apply x-limits via coord_cartesian (does not drop data)
  if (!is.null(limits)) {
    p <- p + coord_cartesian(xlim = limits, expand = FALSE)
  } else {
    p <- p + coord_cartesian(expand = FALSE)
  }

  p <- p + geom_hline(yintercept = 0, linetype = "dashed", color = "grey") 
  p
}


#' check for adequate knots
#' 
#' 
#' 
#' 

check_res_edf <- function(model,k1, k2) {
  res <- residuals(model)
  data <- model.frame(model)
  sm_terms <- gratia::smooths(model)
  te1 <- sm_terms[[1]] |> gsub("te\\(|\\)", "", .)
  te2 <- sm_terms[[2]]
  te3 <- sm_terms[[3]]
  # if k1 and k2 are not provided, extract from model
  if (missing(k1) | missing(k2)) {
    k1 <- model$smooth[[1]]$margin[[1]]$bs.dim
    k2 <- model$smooth[[1]]$margin[[2]]$bs.dim
      } else {
        k1 <- k1
        k2 <- k2
      }
  # construct formula using the actual column names
  form_str <- paste0("res ~", te1,  "),k = c(", 
  k1,",", k2, "), bs = c('tp','cr')) + te(", te2," ),k = c(", 
  k1,",", k2, "), bs = c('tp','cr')) + te(", te3, "), by =weights_lsj, k = c(", 
  k1,",", k2, "), bs = c('tp','cr'))")

  gmod <- mgcv::gam(as.formula(form_str), data = data, gamma = 1.4)
    gmod %>%
    gratia::edf() %>% 
    dplyr::mutate(edf = round(edf, 1))
}
check_res_edf(m_vpd_median)

res <- residuals(m_vpd_median)
gam(res ~ te(vpd, L, by= weights_btl, k = c(3,4), bs = c("tp","cr")), 
 gamma = 1.4, data = moisture_data  |>  select(-c(mean_ndvi, mean_ndmi)) |> drop_na()) |> gratia::edf()


sm_terms <- gratia::smooths(m_vpd_median)

term_variables(m_vpd_median)

sm_terms[[1]] <- {
  s <- sm_terms[[1]]
  if (!grepl("^te\\(", s)) {
    colon <- regexpr(":", s)
    if (colon[1] > 0) {
      inside <- substr(s, 3, colon[1] - 1)   # after "te" up to before ":"
      rest <- substr(s, colon[1], nchar(s))
      s <- paste0("te(", inside, ")", rest)
    } else {
      s <- paste0("te(", substr(s, 3, nchar(s)), ")")
    }
  } else {
    # ensure there is a closing ")" before any ":" (fix malformed cases)
    colon <- regexpr(":", s)
    if (colon[1] > 0 && substr(s, colon[1] - 1, colon[1] - 1) != ")") {
      inside <- substr(s, 4, colon[1] - 1)   # inside of existing "te("
      rest <- substr(s, colon[1], nchar(s))
      s <- paste0("te(", inside, ")", rest)
    }
  }
  s
}


## Simulate some data ....
library(mgcv)
set.seed(1) 
dat <- gamSim(1,n=400,scale=2)

## fit a GAM with quite low `k'
b<-gam(y~s(x0,k=6)+s(x1,k=6)+s(x2,k=6)+s(x3,k=6),data=dat)
plot(b,pages=1,residuals=TRUE) ## hint of a problem in s(x2)

## the following suggests a problem with s(x2)
gam.check(b)

## Another approach (see below for more obvious method)....
## check for residual pattern, removeable by increasing `k'
## typically `k', below, chould be substantially larger than 
## the original, `k' but certainly less than n/2.
## Note use of cheap "cs" shrinkage smoothers, and gamma=1.4
## to reduce chance of overfitting...
rsd <- residuals(b)
gam(rsd~s(x0,k=40,bs="cs"),gamma=1.4,data=dat) ## fine
gam(rsd~s(x1,k=40,bs="cs"),gamma=1.4,data=dat) |>gratia::edf() ## fine
gam(rsd~s(x2,k=10,bs="cs"),gamma=1.4,data=dat)$smooth$bs.dim
gam(rsd~s(x3,k=40,bs="cs"),gamma=1.4,data=dat)  ## fine



get_k_per_smooth <- function(mod) {
  # mod: a fitted mgcv::gam or mgcv::bam object
  out <- lapply(seq_along(mod$smooth), function(j) {
    sm <- mod$smooth[[j]]
    nm <- sm$term             # character vector of covariate names for this smooth
    lab <- sm$label           # nice label
    if (!is.null(sm$margin)) {
      # Tensor product smooth (te/ti/t2): k per margin + product
      km <- vapply(sm$margin, function(mg) mg$bs.dim, integer(1))
      data.frame(
        smooth = lab,
        type   = sm$class,          # e.g. "te", "ti", "t2"
        vars   = paste(nm, collapse = ","),
        k_margins = paste(km, collapse = " x "),
        k_product = prod(km),
        stringsAsFactors = FALSE
      )
    } else {
      # 1D s() smooth
      data.frame(
        smooth = lab,
        type   = sm$class,          # e.g. "s"
        vars   = paste(nm, collapse = ","),
        k_margins = sm$bs.dim,
        k_product = sm$bs.dim,
        stringsAsFactors = FALSE
      )
    }
  })
  do.call(rbind, out)
}

get_k_per_smooth(m_vpd_median)
m_vpd_median$smooth[[1]]$margin[[1]]$bs.dim
m_vpd_median$smooth[[1]]$margin[[2]]$bs.dim
