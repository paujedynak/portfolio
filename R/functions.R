#' ANCOVA model
#'
#' @param batch_factor A string defining the name of the categorical variable
#' @param data A data.frame with input data in a wide format
#'
#' @return
#' @export
batch_effect_ancova <- function(batch_factor, 
                                data) {
  
  formula <- as.formula(paste0("value ~ ", batch_factor))
  aov_fit <- stats::aov(formula, data = data)
  ancova_fit <- car::Anova(aov_fit, type = "III")
  
  return(ancova_fit)
}

#' Compute ANCOVA overall effect
#'
#' @param model_data_long A data frame in a long format
#' @param b_factor A string defining the name of the categorical variable
#'
#' @return
#' @export
compute_ancova <- function(model_data_long, b_factor) {
  
  ancova_tidy <- model_data_long %>%
    group_by(compound_name) %>%
    group_modify(~broom::tidy(batch_effect_ancova(batch_factor = b_factor, .))) %>% 
    ungroup() %>% 
    filter(term == b_factor) %>% 
    dplyr::mutate(p_val_main_eff = p.value,
                  stat_main_eff = statistic) %>% 
    mutate_if(is.numeric, round, 2) %>%
    select(compound_name, term, stat_main_eff, p_val_main_eff)
  
  return(ancova_tidy)
}


#' Linear regression model
#'
#' @param batch_factor A string defining the name of the categorical variable
#' @param data A data.frame with input data in a wide format
#'
#' @return
#' @export
batch_effect_lm <- function(batch_factor, 
                            data) {
  
  formula <- as.formula(paste0("value ~ ", batch_factor))
  lmod <- lm(formula, data = data)
  
  return(lmod)
}


#' Compute linear regression estimates
#'
#' @param model_data_long A data frame in a long format
#' @param b_factor A string defining the name of the categorical variable
#'
#' @return
#' @export
compute_lm <- function(model_data_long, b_factor) {
  
  lmod_tidy <- model_data_long %>%
    group_by(compound_name) %>%
    group_modify(~broom::tidy(batch_effect_lm(batch_factor = b_factor, .), conf.int = TRUE)) %>% 
    ungroup() %>%
    filter(term != "(Intercept)") %>% 
    mutate_if(is.numeric, round, 2) %>%
    rowwise() %>% 
    mutate(est_ci = str_c(estimate, " (", conf.low, "; ", conf.high, ")"),
           p_val = p.value) %>% 
    ungroup() %>% 
    select(compound_name, term, est_ci, p_val)
  
  return(lmod_tidy)
}

#' Correlation equation
#'
#' @param batch_factor 
#' @param data 
#'
#' @return
#' @export
batch_effect_cor <- function(batch_factor, 
                             data) {
  
  corr_fit <- cor(data[["value"]], 
                  data[[batch_factor]],
                  use = "complete.obs", 
                  method = "spearman")
  
  return(corr_fit)
}

#' Compute correlations
#'
#' @param model_data_long A data frame in a long format
#' @param b_factor A string defining the name of the categorical variable
#'
#' @return
#' @export
compute_corr <- function(model_data_long, b_factor) {
  
  corr_effect_pvalue <- model_data_long %>%
    group_by(compound_name) %>%
    group_modify(~skimr::skim(batch_effect_cor(batch_factor = b_factor, .))) %>% 
    ungroup() %>% 
    mutate(corr = numeric.mean) %>% 
    select(compound_name, corr)
  
  return(corr_effect_pvalue)
}


#' Function to fill in data between two values using (Helsel 1990) fill-in method
#'
#' @param to_fill vector of values to fill, should have normal distribution
#' @param lim_inf vector lower bound values, same length as var_to_fill
#' @param lim_sup vector upper bound values, same length as var_to_fill
#' @param seed a scalar defining a seed
#'
#' @return
#' @export
fill_in <- function(to_fill, lim_inf, lim_sup, seed){
  
  # compute percent data below LOD
  pct_below_lod <- sum(to_fill < lim_sup, na.rm = TRUE) / sum(!is.na(to_fill))
  # print(sum(to_fill < lim_sup, na.rm = TRUE))
  
  # warning if pct censored > 30%
  if (pct_below_lod > 0.3) {
    pct <- round(pct_below_lod, 2) * 100
    print(str_c(paste(group, collapse = " - "), ": ", pct,
                "% samples below LOD, should consider categorizing"))
  }
  
  # flag values to impute
  to_impute <- to_fill >= lim_inf & to_fill < lim_sup
  to_impute <- ifelse(is.na(to_impute), as.logical("FALSE"), to_impute)
  
  # flag censored values
  censored <- to_fill < lim_sup 
  
  # compute distribution parameters
  stats_ros <- NADA::cenros(to_fill, censored, forwardT = NULL) #"forwardT = NULL": no log transform
  dist_mean <- NADA::mean(stats_ros)
  dist_sd   <- NADA::sd(stats_ros)
  
  # compute fill-in values by doing a random sample between 0 and LOD for a
  # normal distribution with previsouly computed parmameters
  
  set.seed(seed)
  
  fill_in_vals <- msm::rtnorm(
    n     = sum(to_impute), 
    mean  = dist_mean, 
    sd    = dist_sd, 
    lower = lim_inf, 
    upper = lim_sup
  )
  
  # replace values below LOD by fill in values
  filled_in_var <- to_fill
  filled_in_var[to_impute] <- fill_in_vals
  
  return(filled_in_var)
}

#' A wrapper for the fill_in function when both LOD and LOQ provided
#'
#' @param to_fill vector of values to fill, should have normal distribution
#' @param lim_LOD vector LOD values, same length as var_to_fill
#' @param lim_infer vector lower bound values, same length as var_to_fill
#' @param lim_super vector upper bound values, same length as var_to_fill
#' @param seed a scalar defining a seed
#'
#' @return
#' @export
#'
#' @examples
fun_imp <- function(to_fill, lim_LOD, lim_infer, lim_super, seed) {
  
  imp_to_fill <- fill_in(to_fill = to_fill, lim_inf = lim_LOD, lim_sup = lim_infer, seed = seed)
  imp_to_fill <- fill_in(to_fill = imp_to_fill, lim_inf = lim_infer, lim_sup = lim_super, seed = seed)
  
  return(imp_to_fill)
  
}
