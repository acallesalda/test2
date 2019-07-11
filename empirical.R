## Computing Empirical Power and size examples

# In this file we compute empirical power of the test when the first
# sample is generated from the model 0 and the second sample is
# generated from the model 3.

source("generate_curves.R")
source("DD-plot-test.R")

tests <- logical(length = 100)
for (i in 1:100){
  print('Iteration')
  print(i)
  S <- GenerateCurves()
  S0 <- S[[1]]
  S1 <- S[[4]]
  S0_2 <- S[[1]]
  J <- fdata(S0)
  G <- fdata(S1)
  J_2 <- fdata(S0_2)
  res_pow <- Tester(J, G, B=1000, depth.function = depth.mode, nc = 8)
  res_size <- Tester(J, J_2, B=1000, depth.function = depth.mode, nc = 8)
}

empirical_power <- 1 - sum(res_pow)/100
empirical_size <- 1 - sum(res_size)/100

print('Power')
print(empirical_power)
print('Size')
print(empirical_size)