rm(list = ls())
graphics.off()
cat('\014')
library(dlm)
# data.csv is the time domain PVP signal
pvp = scan(file = 'data.csv', sep = ',')
# initializing DLM
fn <- function(parm) {
  dlmModPoly(order = 1, dV = 0, dW = exp(parm))
}
# finding out MLE using L-BFGS-B
fit <- dlmMLE(pvp, 2, build = fn, method = "L-BFGS-B")
# make sure if converged
(conv <- fit$convergence)
# 
mod <- fn(fit$par)

(FF(mod))
(V(mod))
(GG(mod))
(W(mod))
(m0(mod))
(C0(mod))


