## Usage Examples

source("generate_curves.R")
source("DD-plot-test.R")
source("homogeinity_flores2018.R")
S <- GenerateCurves()
S0 <- S[[1]]
S1 <- S[[4]]
J <- fdata(S0)
G <- fdata(S1)
res_flores <- Tester_Flores(J, G, B=1000, stat=P2, depth.function = depth.mode, nc=8)
res_dd <- Tester(J, G, B=1000, depth.function = depth.mode, nc = 8)