# set of functions and themes commonly used in my plots ####

# paths/directories ####
source('code/00_paths.R')
dir_figs <- paste0(dir_repo, 'final_results/figures/')
dir_tbls <- paste0(dir_repo, 'final_results/tables/')
dir.create(dir_figs, showWarnings = FALSE, recursive=TRUE)
dir.create(dir_tbls, showWarnings = FALSE, recursive=TRUE)

CI95_z <- qnorm(1 - (1 - 0.95)/2)
dash_color <- 'gray50'
# themes ####

theme_set(theme_bw()) # default theme is yuck

# shared theme for multiple tile plots
tileplot_theme <- theme(
  legend.key.size = unit(2, "mm"),
  legend.title = element_text(size=6),
  legend.text = element_text(size=5),
  legend.margin = margin(0,-1,0,-2, unit="mm"),
  legend.text.align = 1,
  plot.title = element_text(size=7),
  plot.subtitle = element_text(size=5),
  axis.title = element_text(size=6),
  axis.text = element_text(size=3.5 , lineheight=0.8)
)

# functions ####

# formats a correlation object as: r = 0.120 (p = 3.0x10^-7)
# using bquote() (set parse=FALSE, such as inside annotate)
annotate_cor.test <- function(
    cor_obj,
    r_digits = 3,
    p_digits = 2,
    p_limit1 = 0.001,
    p_limit2 = 2.2e-16
) {
  r <- cor_obj$estimate
  p <- cor_obj$p.value
  
  r_expr <- get_r_expr(r, r_digits)
  p_expr <- get_p_expr(p, p_digits, p_limit1, p_limit2)
  
  return( bquote(.(r_expr) ~ "(" * .(p_expr) * ")") )
}

# formats a correlation value as: r = 0.120
# using bquote() (set parse=FALSE, such as inside annotate)
get_r_expr <- function(r, r_digits = 3) {
  r_val <- formatC(r, digits = r_digits, format = "f")
  r_expr <- bquote(r == .(r_val))
  return(r_expr)
}

# formats a p-value as: p = 3.0x10^-7
# using bquote() (set parse=FALSE, such as inside annotate)
get_p_expr <- function(p, p_digits = 2, p_limit1 = 0.001, p_limit2 = 2.2e-16) {
  comp <- '=='
  if (p > p_limit1) {
    p_val <- formatC(p, digits = p_digits, format = "fg", flag='#')
    p_expr <- bquote(p == .(p_val))
  } else {
    if (p == 0) {
      p <- p_limit2
      comp <- "<"
    }
    p_exp <- floor(log10(p))
    p_base <- p / 10^p_exp
    p_val <- formatC(p_base, digits = p_digits - 1, format = "f")
    if (comp == '==') {p_expr <- bquote(p == .(p_val) %*% 10^.(p_exp))
    } else {p_expr <- bquote(p < .(p_val) %*% 10^.(p_exp))}
    
  }
  
  return(p_expr)
}
