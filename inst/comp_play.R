library(pika)
library(socorro)
library(viridis)


foo <- function(x) {
    -1 + 2 / (1 + exp(-1 * x))
}

curve(foo(x), from = 0, to = 100)

deathP <- function(x, s, d) {
    dij <- as.matrix(dist(x))
    d <- foo(s)
    # d <- 0
    q <- d + (1 - d) * rowMeans(exp(-(dij / s)^2))

    return(q / sum(q))
}

nspp <- 50

set.seed(30)
x <- rep(sort(rnorm(nspp, 0, 1)), rtnegb(nspp, 1, 0.1))

ss <- exp(seq(log(0.1), log(10), length.out = 5))

pp <- lapply(ss, function(s) deathP(x, s = s, d = 0))

col <- viridis(length(ss))


layout(matrix(c(1, 2, 3, 3), nrow = 2), widths = c(3, 1))
par(oma = c(3.5, 0, 0, 0), mar = c(0, 3, 0, 0) + 0.5, mgp = c(1.75, 0.35, 0), tcl = -0.25)

plot(table(x), xlim = range(x), xaxt = 'n', ylab = 'Abundance', col = 'gray')
abline(v = mean(x), lwd = 3, lty = 3)
text(mean(x), 30, labels = 'mean trait value', pos = 4)

plot(x, pp[[1]], col = col[1], type = 'l', lwd = 2, ylim = range(unlist(pp)),
     xlab = 'Trait values', ylab = 'Death probability')

for(i in 2:length(pp)) {
    lines(x, pp[[i]], col = col[i], lwd = 2)
}

abline(v = mean(x), h = 1 / length(x), lwd = 3, lty = 3)
text(-2.5, 1 / length(x), labels = 'neutral', pos = 3)

par(mar = c(6, 1.5, 6, 3))
plot(rep(1, length(ss)), ss, log = 'y', axes = FALSE, frame.plot = TRUE,
     type = 'n', xlab = '')
segments(x0 = 0.8, x1 = 1.2, y0 = ss, lwd = 3, col = col)
logAxis(4)
mtext(expression(sigma), side = 4, las = 2, cex = 1.5, line = 1.5)
