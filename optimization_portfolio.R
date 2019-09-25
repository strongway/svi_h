library(quantmod)
myStocks <- c("AAPL", "XOM", "GOOGL", "MSFT", "GE", "JNJ", "WMT", "CVX", "PG", "WFC")
getSymbols(myStocks, src = "yahoo")
returns <- lapply(myStocks, function(s)
monthlyReturn(eval(parse(text = s)),
                subset = "2016::2017"))
returns <- do.call(cbind,returns)
colnames(returns) <- myStocks

 library(timeSeries)
 plot(as.timeSeries(returns), at = "chic", minor.ticks="month",
       mar.multi = c(0.2, 5.1, 0.2, 1.1), oma.multi = c(4, 0, 4, 0),
       col = .colorwheelPalette(10), cex.lab = 0.8, cex.axis = 0.8)
 title("Portfolio Returns")
 
  nStocks <- ncol(returns) # number of portfolio assets
  R <- colMeans(returns) # average monthly returns
  S <- cov(returns) # covariance matrix of monthly returns
  s <- sqrt(diag(S)) # volatility of monthly returns
  plot(s, R, type = "n", panel.first = grid(),
        xlab = "Std. dev. monthly returns", ylab = "Average monthly returns")
  text(s, R, names(R), col = .colorwheelPalette(10), font = 2)
  weights <- function(w) # normalised weights
 { drop(w/sum(w)) }
  ExpReturn <- function(w) # expected return
 { sum(weights(w)*R) }
  VarPortfolio <- function(w) # objective function
 {
   w <- weights(w)
   drop(w %*% S %*% w)
 }
 
  fitness <- function(w) # fitness function
 {
   ER <- ExpReturn(w)-0.01
   penalty <- if(ER < 0) 100*ER^2 else 0
   -(VarPortfolio(w) + penalty)
 }
 
    GA <- ga(type = "real-valued", fitness = fitness,
              min = rep(0, nStocks), max = rep(1, nStocks), names = myStocks,
              maxiter = 1000, run = 200, optim = TRUE)
  summary(GA)
 (w <- weights(GA@solution))
 
  ExpReturn(w)

  VarPortfolio(w)

  
  
  barplot(w, xlab = "Stocks", ylab = "Portfolio weights",
           cex.names = 0.7, col = .colorwheelPalette(10))
 
  