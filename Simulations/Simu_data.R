#-------------Simulation de jeux de données-------------------

local.dir = getwd()
sourceCpp(file.path(local.dir,"src/optim_rank_ZIP.cpp"))

n <- 300
p <- 30
d <- 2
q <- 3

X <- cbind(c(rep(1, n*p)), matrix(rnorm(n*p*d), nrow = n*p))

B <- c(0.5, rnorm(d))
B
D <-  c(2, -1, 3)

mu <- VectorToMatrix(X%*%B, n, p)
nu <- VectorToMatrix(X%*%D, n, p)


Prob <- plogis(nu)
U <- matrix(rbinom(n*p, p = Prob, size = 1), nrow = n)

plot(nu, U)
lines(sort(nu), plogis(sort(nu)))
boxplot(nu ~ U)

W <- matrix(rnorm(n*q), nrow = n)
C <- matrix(rnorm(p*q,0,1), nrow = p)/ sqrt(q)

var(X%*%B)
var(MatrixToVector(W%*%t(C)))

Lambda <- exp(mu + W %*% t(C))
Z <- matrix(rpois(n*p, lambda = Lambda), nrow = n)

Y <- ifelse(U == 0, 0, Z)
range(Y)

Y.na10 <- prodNA(Y, 0.1)
Y.na30 <- prodNA(Y.na10, 0.3)
Y.na60 <- prodNA(Y.na30, 0.6)

Data <- list(Y = Y, Z = Z, C = C, B = B, X = X, W = W, mu = mu, nu = nu, D = D, Y.na10 = Y.na10, Y.na30 = Y.na30, Y.na60 = Y.na60)
save(Data, file = "Simu_data.Rdata")
