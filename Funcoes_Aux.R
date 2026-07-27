# ******************************************************************************
# SCRIPT DE FUNÇÕES AUXILIARES PARA AVALIAÇÕES IMOBILIÁRAS
# AUTOR: LUIZ FERNANDO PALIN DROUBI (lfpdroubi@gmail.com)
# ******************************************************************************

# ==============================================================================
## MÉDIA HARMÔNICA
# ==============================================================================

hmean <- function(x, w, ...) {
  n <- length(x)
  inv_x <- lapply(x, FUN = function(x) 1/x)
  if (missing(w)) {
    hm <- n/Reduce("+", inv_x)
  } else {
    inv_wx <- mapply(FUN = `*`, inv_x, w)
    hm <- sum(w)/rowSums(inv_wx)
  }
  return(hm)
}

# ==============================================================================
## COEFICIENTE DE DETERMINAÇÃO NORMAL E AJUSTADO
# ==============================================================================

R2 <- function(model) {
  y <- model$fitted.values + model$residuals
  res <- model$residuals
  w <- model$weights
  
  if (is.null(w)) {
    TSS <- sum((y - mean(y))^2)
    RSS <- sum(res^2)
  } else {
    TSS <- sum(w*(y-weighted.mean(y,w))^2)
    RSS <- sum(w*(res)^2)
  }
  
  R2 <- 1 - RSS/TSS
  return(R2)
}

adjR2 <- function(model) {
  R2 <- R2(model)
  n <- nobs(model)
  gl <- df.residual(model)
  
  R2adj <- adjR2.default(R2 = R2, n = n, dof = gl)
  return(R2adj)
}

adjR2.default <- function(R2, n, dof) {
  R2adj <- 1-(1-R2)*(n-1)/dof
  return(R2adj)
}

function(model) {
  res <- residuals(model)
  hii <- hatvalues(model)
  y <- model$model[, 1]
  w <- model$weights
  pr <- res/(1-hii)
  
  if (is.null(w)) {
    TSS <- sum((y - mean(y))^2)
    PRESS <- sum(pr^2)
  } else {
    TSS <- sum(w*(y-weighted.mean(y,w))^2)
    PRESS <- sum(w*(pr)^2)
  }
  
  pred.r.sqr <- 1 - PRESS/TSS
  return(pred.r.sqr)
}

# ==============================================================================
## FUNÇÕES PARA FORMATAÇÃO DE TEXTO
# ==============================================================================

brf <- function(x, decimal.mark = ",", big.mark = ".", digits = 2,
                nsmall = 2, scientific = FALSE, ...) {
  format(x, decimal.mark = decimal.mark, big.mark = big.mark, digits = digits,
         nsmall = nsmall, scientific = scientific, ...)
}


pct <- function (x, ...) {
  if (length(x) == 0)
    return(character())
  x <- brf(100*x, ...)
  paste0(x, "%")
}

# ==============================================================================
## FUNÇÕES ESPACIAIS
# ==============================================================================

cria_grid <- function(limite_polygon, cell_size = 20)
{
  gridPpts <- sf::st_make_grid(limite_polygon,
                               what = "centers",
                               cellsize = c(cell_size, cell_size))
  gridPpts <- sf::st_sf(gridPpts)
  
  return(gridPpts[, ncol(gridPpts)])
}