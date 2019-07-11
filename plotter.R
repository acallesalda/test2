#N = 50

pot1 <- c(0.05, 0.05, 0.2, 0.29, 0.67, 0.88, 0.93, 0.99, 1) #pot para n = 50, mi test
deltas <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4) #pot para n = 50, flores
pot2 <- c(0.07, 0.08, 0.25, 0.23, 0.52, 0.67, 0.81, 0.88, 0.90)
plot(deltas, pot1, 'o', pch=1, xlab=expression(delta), ylab='Empircal Power', main='Empirical Power for both tests with n = 50')
lines(deltas, pot2, 'o', pch=2)
legend("topleft", legend=c('DD-plot', 'Flores (2018)'), pch=c(1,2))

#mierda

# N = 250

pot1 <- c(0.05, 0.26, 0.68, 0.98, 1)
pot2 <- c(0.05, 0.15, 0.35, 0.65, 0.80)
#falta pot2
deltas <- c(0, 0.05, 0.1, 0.15, 0.2)
plot(deltas, pot1, 'o', pch=1, xlab=expression(delta), ylab='Empircal Power', main='Empirical Power for both tests with n = 50')
lines(deltas, pot2, 'o', pch=2)
legend("topleft", legend=c('DD-plot', 'Flores (2018)'), pch=c(1,2))


# n = 50
#Suficiente de mu, ahora necesito cambiar los de la varianza

pot1 <- c(0.05, 0.44, 0.50, 0.65, 0.67, 0.83, 0.90, 0.90, 0.91, 0.94, 0.92, 0.94, 0.97, 0.94, 0.99, 1)
pot2 <- c(0.05, 0.13, 0.16, 0.25, 0.47, 0.40, 0.51, 0.56, 0.61, 0.68, 0.65, 0.66, 0.58, 0.70, 0.65, 0.69)
cons <- c(0.3, 0.9, 1.8,  2.7,  3.6,  4.5,  5.4,  6.3,  7.2,  8.1,  9.0,  9.9, 10.8, 11.7, 12.6, 13.5)
plot(cons, pot1, 'o', pch=2, xlab='k', ylab='Empircal Power', main='Empirical Power for both tests with n = 50')
lines(cons, pot2, 'o', pch=1)
legend("topleft", legend=c('Flores (2018)','DD-plot'), pch=c(1,2), cex = 0.8)

#N = 20
pot1 <- c(0.05, 0.27, 0.34, 0.43, 0.32, 0.53, 0.52, 0.56,0.51, 0.55, 0.53, 0.59, 0.55, 0.72, 0.66, 0.70, 0.66, 0.64, 0.70, 0.74, 0.74)
cons <- c(0.3, 1.5,  3.0,  4.5,  6.0,  7.5,  9.0, 10.5, 12.0, 13.5, 15.0, 16.5, 18.0, 19.5, 21.0, 22.5, 24.0, 25.5, 27.0, 28.5, 30.0)
pot2 <- c(0.05, 0.20, 0.32, 0.37, 0.47, 0.40, 0.52, 0.51, 0.61, 0.57, 0.60, 0.59, 0.61, 0.61, 0.64, 0.61, 0.64, 0.67, 0.62, 0.69, 0.66)
plot(cons, pot1, 'o', pch=2, xlab='k', ylab='Empircal Power', main='Empirical Power for both tests with n = 20')
lines(cons, pot2, 'o', pch=1)
legend("topleft", legend=c('Flores (2018)','DD-plot'), pch=c(1,2), cex = 0.75)


#N = 250
pot1 <- c(0.05, 0.12, 0.41, 0.61, 0.78, 0.89, 0.93, 0.91, 0.97, 0.99, 0.97)
cons <- c(0.3, 0.345, 0.390, 0.435, 0.480, 0.525, 0.570, 0.615, 0.660, 0.705, 0.750)
pot2 <- c(0.05, 0.10, 0.06, 0.06, 0.07, 0.09, 0.06, 0.14, 0.12, 0.13, 0.10)
plot(cons, pot1, 'o', pch=2, xlab='k', ylab='Empircal Power', main='Empirical Power for both tests with n = 250')
lines(cons, pot2, 'o', pch=1)
legend("topleft", legend=c('Flores (2018)','DD-plot'), pch=c(1,2), cex = 0.75)
