library(magrittr)

make_me_a_biplot <- function(data.row, data.col){
  df <- data.frame(x = c(data.row[,1], data.col[,1]),
                   y = c(data.row[,2], data.col[,2]),
                   z = c(paste(rep("Row", nrow(data.row)), as.character(1:nrow(data.row)), sep = ""),
                         paste(rep("Col", nrow(data.col)), as.character(1:nrow(data.col)), sep = "")))
  margins <- par()$mar
  par(mar = c(3, 3, 1, 1))
  plot(1, type = "n", xlab = "", ylab = "")
  text(df$x, df$y, labels = df$z)
  segments(x0 = rep(0, nrow(data.row)), y0 = rep(0, nrow(data.row)), 
           x1 = df$x[1:nrow(data.row)], y1 = df$y[1:nrow(data.row)], col = "red")
  segments(x0 = rep(0, nrow(data.col)), y0 = rep(0, nrow(data.col)), 
           x1 = df$x[nrow(data.row)+1:nrow(data.row)+nrow(data.col)], 
           y1 = df$y[nrow(data.row)+1:nrow(data.row)+nrow(data.col)], lty = 2, col = "red")
  abline(h = 0, lty = 3)
  abline(v = 0, lty = 3)
  par(mar = margins)
}

#The rank 2 case
(y <- c(rep(2, 2), -4, 2, 1, -3, 0, -3/2, 3/2, -1, -1/2, 3/2) %>% 
  matrix(data = ., nrow = 4, byrow = TRUE))
y_svd <- y %>% svd(., nu = 2, nv = 2)
(g <- y_svd$u %*% diag(y_svd$d[1:2]))
(h <- y_svd$v)
df <-  data.frame(x = c(g[,1], h[,1]), 
                  y = c(g[,2], h[,2]), 
                  z = c(paste(rep("Row", 4), as.character(1:4), sep = ""),
                        paste(rep("Col", 3), as.character(1:3), sep = "")))
margins <- par()$mar
par(mar = c(3, 3, 1, 1))
plot(1, type = "n", xlab = "", ylab = "", xlim = c(-5, 2), ylim = c(-1.5, 1))
text(df$x, df$y, labels = df$z)
segments(x0 = rep(0, 4), y0 = rep(0, 4), x1 = df$x[1:4], y1 = df$y[1:4], col = "red")
segments(x0 = rep(0, 3), y0 = rep(0, 3), x1 = df$x[5:7], y1 = df$y[5:7], lty = 2, col = "red")
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)
par(mar = margins)

#A different factorization
(g2 <- y_svd$u)
(h2 <- y_svd$v %*% diag(y_svd$d[1:2]))
df <-  data.frame(x = c(g2[,1], h2[,1]), 
                  y = c(g2[,2], h2[,2]), 
                  z = c(paste(rep("Row", 4), as.character(1:4), sep = ""),
                        paste(rep("Col", 3), as.character(1:3), sep = "")))
margins <- par()$mar
par(mar = c(3, 3, 1, 1))
plot(1, type = "n", xlab = "", ylab = "", xlim = c(-3.5, 5.5), ylim = c(-1.3, 1))
text(df$x, df$y, labels = df$z)
segments(x0 = rep(0, 4), y0 = rep(0, 4), x1 = df$x[1:4], y1 = df$y[1:4], col = "red")
segments(x0 = rep(0, 3), y0 = rep(0, 3), x1 = df$x[5:7], y1 = df$y[5:7], lty = 2, col = "red")
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)
par(mar = margins)

#Scaling
y_svd2 <- y %>% scale(., center = FALSE, scale = TRUE) %>% svd(., nu = 2, nv = 2)
(g3 <- y_svd2$u %*% diag(y_svd2$d[1:2]))
(h3 <- y_svd2$v)
df <-  data.frame(x = c(g3[,1], h3[,1]), 
                  y = c(g3[,2], h3[,2]), 
                  z = c(paste(rep("Row", 4), as.character(1:4), sep = ""),
                        paste(rep("Col", 3), as.character(1:3), sep = "")))
margins <- par()$mar
par(mar = c(3, 3, 1, 1))
plot(1, type = "n", xlab = "", ylab = "", xlim = c(-3, 1), ylim = c(-1, 1))
text(df$x, df$y, labels = df$z)
segments(x0 = rep(0, 4), y0 = rep(0, 4), x1 = df$x[1:4], y1 = df$y[1:4], col = "red")
segments(x0 = rep(0, 3), y0 = rep(0, 3), x1 = df$x[5:7], y1 = df$y[5:7], lty = 2, col = "red")
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)
par(mar = margins)

#The rank R case
(y_big <- c(1:4, 1, 0, -1, 2, 4, 2, rep(0, 2), 4, rep(0, 2), seq(-3, -12, -3), 4, rep(0, 2), rep(1, 3)) %>% 
    matrix(data = ., nrow = 5, byrow = TRUE))
qr(y_big)$rank
y_big_svd <- y_big %>% svd(., nu = 2, nv = 2)
#Goodness of fit
sum(y_big_svd$d[1:2]^2)/sum(y_big_svd$d^2)
(g <- y_big_svd$u %*% diag(y_big_svd$d[1:2]))
(h <- y_big_svd$v)
df <-  data.frame(x = c(g[,1], h[,1]), 
                  y = c(g[,2], h[,2]), 
                  z = c(paste(rep("Row", 5), as.character(1:5), sep = ""),
                        paste(rep("Col", 5), as.character(1:5), sep = "")))
margins <- par()$mar
par(mar = c(3, 3, 1, 1))
plot(1, type = "n", xlab = "", ylab = "", xlim = c(-6, 20), ylim = c(-4, 1))
text(df$x, df$y, labels = df$z)
segments(x0 = rep(0, 5), y0 = rep(0, 5), x1 = df$x[1:5], y1 = df$y[1:5], col = "red")
segments(x0 = rep(0, 5), y0 = rep(0, 5), x1 = df$x[6:10], y1 = df$y[6:10], lty = 2, col = "red")
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)
par(mar = margins)

#Scale
(y_big <- c(1:4, 1, 0, -1, 2, 4, 2, rep(0, 2), 4, rep(0, 2), seq(-3, -12, -3), 4, rep(0, 2), rep(1, 3)) %>% 
    matrix(data = ., nrow = 5, byrow = TRUE))
y_big_svd2 <- y_big %>% scale(., center = FALSE, scale = TRUE) %>%  svd(., nu = 2, nv = 2)
#Goodness of fit
sum(y_big_svd2$d[1:2]^2)/sum(y_big_svd2$d^2)
(g2 <- y_big_svd2$u %*% diag(y_big_svd2$d[1:2]))
(h2 <- y_big_svd2$v)
df <-  data.frame(x = c(g2[,1], h2[,1]), 
                  y = c(g2[,2], h2[,2]), 
                  z = c(paste(rep("Row", 5), as.character(1:5), sep = ""),
                        paste(rep("Col", 5), as.character(1:5), sep = "")))
margins <- par()$mar
par(mar = c(3, 3, 1, 1))
plot(1, type = "n", xlab = "", ylab = "", xlim = c(-1.5, 4.5), ylim = c(-1.5, 0.5))
text(df$x, df$y, labels = df$z)
segments(x0 = rep(0, 5), y0 = rep(0, 5), x1 = df$x[1:5], y1 = df$y[1:5], col = "blue")
segments(x0 = rep(0, 5), y0 = rep(0, 5), x1 = df$x[6:10], y1 = df$y[6:10], lty = 2, col = "red")
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)
par(mar = margins)

#PCA
data("mtcars")
biplot(prcomp(mtcars, scale = TRUE))

x<-prcomp(mtcars, scale=TRUE)
choices = 1L:2L
scale = 1
pc.biplot = FALSE
scores<-x$x
lam <- x$sdev[choices]
n <- NROW(scores)
lam <- lam * sqrt(n)
lam <- lam^scale
(yy<-t(t(x$rotation[, choices]) * lam))
(xx<-t(t(scores[, choices])/lam))
biplot(xx,yy)
