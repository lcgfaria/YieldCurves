# Identification ----------------------------------------------------------
# 
# FGV - Fundacao Getulio Vargas
# EESP - Escola de Economia São Paulo
# MPEF - Mestrado Profissional em Economia e Finanças
#
# YIELD CURVES COMPARISON
# Autor: Luis Giovanni Faria
#
# Date: 13/09/2022
#

## Functions Library

ReduceSparseMatrix <- function(sparse) {
  ## Description: reduce sparse matrix to only significant columns
  ##   sparse: sparse matrix
  ## Return
  ##   reduced: reduced matrix
  reduced <- data.frame(matrix(nrow=nrow(sparse), ncol=0))
  for (i in 1:ncol(sparse)) {
    if (sum(sparse[,i])!=0) {
      reduced[,ncol(reduced)+1] <- sparse[,i]
      colnames(reduced)[ncol(reduced)] <- i
    } # end-if
  } # end-for 
  return(reduced)
} # end-function


RootMeanSquareError <- function(fi,params,t,yield,w=1) {
  ## Compute RMSE for a given yield function at time t
  ## fi: yield function
  ## params: vector with function parameters
  ## t: vector with times
  ## yield: vector with wields
  ## w: vector with weights [1|1:n]:  if (length(w)<length(t)) w = rep(w,length(t))
  return(sqrt(sum(((yield - fi(params,t))^2)*w)))
}


CCSplineLGF <- function(x,y,xstar) {
  ## Description: Constrained Cubic Spline Interpolation
  ## Arguments:
  ##    x: Numeric vector with x values
  ##    y: Numeric vector with y values
  ##    xstar: x value to interpolate y*Valor de x para fazer a interpola??o: "y*" 
  ## Return
  ##    ystar: y interpolated value
  
  #
  # Verifica se xstar está fora do vetor x
  if ( xstar<x[1] || xstar>x[length(x)] ) {
    ystar <- NA           # fora dos limites do vetor
  } else if (xstar==min(x)) {
    ystar <- y[1]         # primeiro elemento do vetor
  } else if (xstar==x[length(x)]) {
    ystar <- y[length(x)] # último elemento do vetor
  } else {
    ## Procura intervalo de interpolação
    ##  Algoritimo golden search para números inteiros
    r <- (sqrt(5)-1)/2    # razão áurea
    a <- 1                # início do intervalo de pesquisa
    b <- length(x)        # fim do intervalo de pesquisa
    h <- b-a              # tamanho do intervalo de pesquisa
    c <- trunc(a+(1-r)*h) # razão áurea inferior, arredonda c para baixo
    d <- trunc(a+r*h+0.5) # razão áurea superior, arredonda d para cima
    ## Pesquisa enquanto tamanho do intervalo for maior que um
    while ( h>1 ) {
      if (x[d]<=xstar) {
        a <-d
      } else if (xstar<=x[c]) {
        b <- c
      } else {
        a <- c
        b <- d
      }
      h <- b-a              # tamanho do intervalo de pesquisa
      c <- trunc(a+(1-r)*h) # arredonda c para baixo
      d <- trunc(a+r*h+0.5) # arredonda d para cima
    } # end-while
    i <- b
    ## Calcula fórmulas de interpolação
    if (i==length(x)) {
      df_x1_ant <- 2/( (x[i]-x[i-1])/(y[i]-y[i-1]) + (x[i-1]-x[i-2])/(y[i-1]-y[i-2]) ) # Eq 7a
      df_x1 <- ( (3*(y[i]-y[i-1]))/(2*(x[i]-x[i-1])) ) - df_x1_ant/2   # Eq 7b
    } else {
      df_x1 <- 2/( (x[i+1]-x[i])/(y[i+1]-y[i]) + (x[i]-x[i-1])/(y[i]-y[i-1]) ) # Eq 7a  
    } # end-if
    if (i==2) {
      df_x0 <- ( (3*(y[i]-y[i-1]))/(2*(x[i]-x[i-1])) ) - df_x1/2   # Eq 7b
    } else {
      df_x1_ant <- 2/( (x[i]-x[i-1])/(y[i]-y[i-1]) + (x[i-1]-x[i-2])/(y[i-1]-y[i-2]) ) # Eq 7a
      df_x0 <- df_x1_ant
    } # end-if
    ddf_x0 <- -(2*(df_x1+2*df_x0)/(x[i]-x[i-1]))+6*(y[i]-y[i-1])/(x[i]-x[i-1])^2  # Eq 8
    ddf_x1 <- 2*(2*df_x1+df_x0)/(x[i]-x[i-1]) - 6*(y[i]-y[i-1])/(x[i]-x[i-1])^2  # Eq 9    
    d <- (ddf_x1-ddf_x0)/6/(x[i]-x[i-1])  # Eq 10
    c <- (x[i]*ddf_x0-x[i-1]*ddf_x1)/2/(x[i]-x[i-1])  # Eq 11
    b <- ((y[i]-y[i-1])-c*(x[i]^2-x[i-1]^2)-d*(x[i]^3-x[i-1]^3))/(x[i]-x[i-1])  # Eq 12
    a <- y[i]-b*x[i]-c*x[i]^2-d*x[i]^3  # Eq 13
    ystar <- a+b*xstar+c*xstar^2+d*xstar^3  # Eq 1
  } # end-case
  ## retorna ystar
  return(ystar)    
} # end-function


RootNewtonRaphson <- function(f, x0,..., epson=1e-4, nmax=1000) {
  ## Description: Function root using Newton-Raphson method
  ## Arguments:
  ##  f: function to find root
  ##  x0: initial value of x to start search
  ##  epson: precision of root value
  ##  nmax: maximum number of interactions before search is aborted
  ## Return
  ##  xn: root 
  x1 <- x0
  fx <- f(x1,...)
  n <- 1
  while (abs(fx)>epson && n<nmax) {
    dx <- ( f(x1+epson,...) - f(x1-epson,...) ) / (2*epson)
    if (dx==0) { 
      x1 <- x1+epson 
    } else {
      x1 <- x1 - fx/dx
    } # end-else
    fx <- f(x1,...)
    n <- n+1
  } # end-while
  if (abs(fx)<epson) {
    return(x1)
  } else {
    return(NaN)
  } # end-if
} # end-function


PricePU <- function(yield,days,faceValue=100,coupon=0.05,daysYear=365,cc=TRUE,ac=FALSE) {
  ## Description: Compute bond price
  ## Arguments:
  ##  yield: YTM
  ##  days: vector days to coupon payment
  ##  faceValue
  ##  coupon rate semiannual
  ##  cc: if use continuous capitalization of not
  ##  br: consider accrued coupon on first payment
  ## Return
  ##  price
  price <- 0
  n <- length(days)-1
  i <- 1
  firstCoupon=1
  if (length(days)==1&&ac) firstCoupon <- coupon*(daysYear/2-days[i]) else firstCoupon <- 1
  while (i <= n) {
    if (i==1&&ac) firstCoupon <- coupon*(daysYear/2-days[i]) else firstCoupon <- 1
    if (cc) {
      price <- price+faceValue*coupon*firstCoupon*exp(-yield*days[i]/daysYear)
    } else {
      price <- price+faceValue*coupon*firstCoupon/((1+yield)^(days[i]/daysYear));
    } # end-if
    i <- i+1
  } # end-while
  if (cc) {
    price <- price+faceValue*(1+coupon*firstCoupon)*exp(-yield*days[i]/daysYear)  
  } else {
    price <- price+faceValue*(1+coupon*firstCoupon)/((1+yield)^(days[length(days)]/daysYear))
  } # end-else
  return(price)
} # end-function


YieldToMaturity <- function(price,x0=0.0,...) {
  ## Description: Calculate YTM root
  DeltaPricePU <- function(yield,price,...) {
    ## Description: Auxiliary function to calculate YTM root
    return(price-PricePU(yield,...))
  } # end-function
  return(RootNewtonRaphson(DeltaPricePU,x0,price,...))
} # end-function


VarBox <- function(x,lower,upper){
  ## Description: Check if variables are within limits
  ## Arguments:
  ##  x: vector with values to fit in the box
  ##  upper: vector with upper limits
  ##  lower: vector with lower limits
  ## Return
  ##  x: vector with values within box limits   
  for(i in 1:length(x)){
    if(x[i]<lower[i]) x[i] <- lower[i]
    if(x[i]>upper[i]) x[i] <- upper[i]
  } # end-for
  return(x)
} # end-function


DifferentialEvolution = function(fx,fi,xG0,...,F=1,CR=0.5,nmaxG=100,lower=-Inf,upper=Inf,hx=NULL,gx=NULL){
  ## Description: Optimization using Differential Evolution Method
  ## Arguments:
  ##  f: function do me optimized for main function fi
  ##  fi: Main function
  ##  xG0: matrix with initial search
  ##  lower,upper: box limits for xG
  ##  F: mutation factor
  ##  CR: cross over factor
  ## Return
  ##  xMin: vector with optimum point
  ##  
  ## Verify restrictions  
  if (!is.null(hx) || !is.null(gx)){
      fobj <- function(fi,x,...) {
      ## Original function
      s <- fx(fi,x,...)
      ## Set precision
      epson <- 10^(trunc(log(abs(s)+1,10))-7)
      ## Equality estrictions
      if (!is.null(hx)){
        for (i in 1:length(hx)) {
          hxi <- hx[[i]](x)
          s <- s + (hxi^2)/epson
        } # end-if
      } # end-if
      ## Inequality restrictions
      if (!is.null(gx)){
        for (i in 1:length(gx)) {
          gxi <- gx[[i]](x)
          # s <- s + ifelse(gxi<=epson,0,exp(gxi)/epson);
          s <- s + ifelse(gxi<=0,0,exp(gxi)/epson)
        } # end-if
      } # end-if
      return(s)
    } # end-function
  } else {
    fobj <- fx
  } # end-else
  ## Initialize current generation
  xG <- xG0
  ## Population size
  NP <- ncol(xG)
  ## Number of genes
  D <- nrow(xG)
  ## Treat input parameters
  if (length(lower)<D) lower <- rep(lower,D)
  if (length(upper)<D) upper <- rep(upper,D)
  ## Generation cost
  costG <- c()
  ## Calculate cost of current generation
  for (i in 1:NP) {
    ## Fit in the box
    xG[,i] <- VarBox(xG[,i],lower,upper)
    ## Calculate cost
    costG[i] <- fobj(fi,xG[,i],...)
  } # end-for
  ## --- Generation oop
  for (g in 1:nmaxG){
    ## Target-vector loop
    for (i in 1:NP){
      ## --- Mutation
      ## Creates mutant individual: mutant-vector: vi
      r <- sample((1:NP)[-c(i)],3)
      ## Mutant-vector
      viG1 <- xG[,r[1]] + F*(xG[,r[2]] - xG[,r[3]]);
      ## --- Crossover
      ## Initialize trial-vector
      uiG1 <- xG[,i]
      ## rnbr_i: choose mutant vector gene for trial vector
      rnbr_i <- sample(1:D,1)
      for (j in 1:D){
        ## For each gene random value to compare with  CR
        randb_j <- runif(1,min=0,max=1)
        ## Set trial-vector
        if (randb_j <= CR || j == rnbr_i) {
          uiG1[j] <- viG1[j]
        } # end-if
      } # end-for
      ## --- Selection
      ## Fits inside box
      uiG1 <- VarBox(uiG1,lower,upper)
      ## Calculate trial vector score
      cost_uiG1 <- fobj(fi,uiG1,...)
      ## Select most evolved individual
      if (cost_uiG1 < costG[i]){
        xG[,i] <- uiG1
        costG[i] <- cost_uiG1
      } # end-if
    } # end-for
  } ## Generation
  ## Choose most evolved individual with lower cost
  id <- which(costG==min(costG))[1]
  xMin <- xG[,id]
  ## Return optimum vector
  return(xMin)
} # end-function


##########################################################################


Gradient <- function(f,fi,x,...) {
  ## Description: Compute Gradient of function
  ## Arguments:
  ##  f: function
  ##  x: vector with coordinates
  ## Return
  ##  g: vector with gradient   
  g = c()
  for (i in 1:length(x)) {
    h <- 10^(trunc(log(abs(x[i])+1,10))-6)
    xplus <- x
    xminus <- x
    xplus[i] <- x[i]+h
    xminus[i] <- x[i]-h
    g[i] <- (f(fi,xplus,...)-f(fi,xminus,...))/(2*h)
  } # end-for
  return(g)
} # end-function


OptimizeSteepestDescent <- function(f, fi, x0, ..., lower=-Inf, upper=Inf, epson = 1e-4, nmax=10000) {
  ## Description: Optimization using Gradient/Steepest Descent Method
  ## Arguments:
  ##  f: function 
  ##  fi: function
  ##  x0: vector with initial search
  ##  lower,upper: box limits for x
  ##  epson: solution precision
  ##  nmax: max number of search interactions
  ## Return
  ##  y: vector with optimum point
  ##
  ## verify parameters
  if (length(lower)<length(x0)) lower = rep(lower,length(x0))
  if (length(upper)<length(x0)) upper = rep(upper,length(x0))
  ## counter
  n <- 0
  ## Initialize x[n-1]
  xn_1 <- x0
  ## Auxiliary variables
  ## differençe to x
  xdiff <- 1+epson
  ## difference to f
  fdiff <- 1+epson
  while ((xdiff>epson || fdiff>epson) && n<nmax) {
    ## calculate gradient
    gxn_1 <- Gradient(f,fi,xn_1,...)
    ## gradient vector
    ugxn_1 <- gxn_1/sqrt(sum(gxn_1*gxn_1));
    ## step size
    r <- (sqrt(5)-1)/2
    ## limits a and b
    a <- 0
    b <- sqrt(sum(gxn_1*gxn_1))
    h <- b - a
    k <- 0
    while(h>epson && k<nmax){
      ## new limits
      c <- a + (1-r)*h
      d <- a + r*h
      ## box
      xc <- VarBox(xn_1 - c*ugxn_1,lower,upper)
      xd <- VarBox(xn_1 - d*ugxn_1,lower,upper)
      ## new interval
      if (f(fi,xc,...)<f(fi,xd,...)){
        b <- d
      } else {
        a <- c
      } # end-if
      h <- b - a
      k <- k + 1
    } # end-while
    if (f(fi,xc,...)<f(fi,xd,...)){
      alphan_1 <- c
    } else {
      alphan_1 <- d
    } # end-else
    ## calculate next point
    xn <- VarBox(xn_1 - alphan_1*ugxn_1,lower,upper)
    ## update differences
    fdiff <- abs(f(fi,xn,...)-f(fi,xn_1,...))
    xdiff <- sqrt(sum((xn_1-xn)*(xn_1-xn)))
    n <- n + 1
    ## new interaction
    xn_1 <- xn
  } # end-while
  return(xn)  
} # end-function


#' CCSpline: Constrained Cubic Spline
#'
#' @param x [1:n] Valores x (variável independente) da série de dados.*Ordem Crescente*
#' @param y [1:n] Valores y (variável dependente) da série de dados
#' @param xstar [1] Valor x* para o qual se deseja fazer a interpolação
#'
#' @return [1] Valor y* da interpolação
#' @export
#'
#' @examples CCSpline(x=c(0,10,30,50),y=c(30,130,150,150),xstar=20)
#' 
CCSpline = function(x,y,xstar){
  
  ## Variáveis auxiliares
  dfxi = 0;     ## f'(xi)
  dfxi_1 = 0;   ## f'(xi_1)
  dfixi_1 = 0;  ## f'i(xi_1)
  dfixi = 0;    ## f'i(xi)
  ddfixi_1 = 0;  ## f''i(xi_1)
  ddfixi = 0;    ## f''i(xi)
  
  ## Coeficientes da equação
  ai = 0;
  bi = 0;
  ci = 0;
  di = 0;
  
  ## Atencao: 
  ## No paper: i = 0,1,2,...n-1
  ## No R: contadorex de x i = 1,2,3,...n, ou seja, o primeiro i=1
  ## No R: contadores de f i = 2, 3, 4,...n. Ponto x2: f2 a esquerda de x2 e f3 a direita de x2
  
  ## Numero de pontos = índice da última função
  n = length(x);
  
  ## Contador das funções
  i = 2;
  
  ## Calcular apenas para os segmentos necessarios, ou seja, xi_1 < xstar < x_i
  while (x[i-1] <= xstar && i <= n){
    
    ## -- Calcula f'(xi) = dfxi
    
    ## Verifica se é o ultimo ponto
    if (i < n) {
      ## Slope a esquerda de xi
      coefLeft = (y[i]-y[i-1])/(x[i]-x[i-1]);
      
      ## Slope a direita de xi
      coefRight = (y[i+1]-y[i])/(x[i+1]-x[i]);
      
      ## Verifica se houve mudanca de sinal do slope e se é diferente de zero
      ## Verificação alternativa ((y[i]-y[i-1]) * (y[i+1]-y[i])) > 0
      if (sign(coefLeft)==sign(coefRight) && coefLeft!=0){
        dfxi = 2/((x[i+1]-x[i])/(y[i+1]-y[i])+(x[i]-x[i-1])/(y[i]-y[i-1]))
      } else {
        ## Se muda o sinal ou coefLeft==coefRight==0
        dfxi = 0;
      }
      
      ## Para i < n temos f'i(xi)=f'(xi)
      dfixi = dfxi;
      
    } else {
      # Para i == n temos f'i(xi_1)=[f'(xi) da iteracao anterior]
      dfxi_1 = dfxi;
    }
    
    ## Se for a primeira funcao i = 2 calcula f'2(x1)
    if (i == 2){
      dfixi_1 = (3*(y[i]-y[i-1])/(2*(x[i]-x[i-1])))-dfxi/2;
    }
    
    ## Se for a ultima funcao i = n calcula f'i(xi)=f'n(xn)
    if (i == n){
      dfixi = (3*(y[i]-y[i-1])/(2*(x[i]-x[i-1])))-dfxi_1/2;
    }
    
    ## Se o próximo x[i]>=xstar então calcula o polinomio
    if (x[i]>=xstar){
      
      ## Calcula das derivadas de segunda ordem
      ddfixi_1 = - 2*(dfixi+2*dfixi_1)/(x[i]-x[i-1])+6*(y[i]-y[i-1])/((x[i]-x[i-1])^2)
      
      ddfixi = 2*(2*dfixi+dfixi_1)/(x[i]-x[i-1])-6*(y[i]-y[i-1])/((x[i]-x[i-1])^2)
      
      ## Calcula os coeficientes
      
      di = (ddfixi - ddfixi_1)/(6*(x[i] - x[i-1]));
      
      ci = (x[i]*ddfixi_1 - x[i-1]*ddfixi)/(2*(x[i] - x[i-1]));
      
      bi = ((y[i] - y[i-1]) - ci*(x[i]^2-x[i-1]^2) - di*(x[i]^3-x[i-1]^3))/(x[i] - x[i-1]);
      
      ai = y[i-1] - bi*x[i-1] - ci*x[i-1]^2 - di*x[i-1]^3;
      
    }
    
    ## Atualiza para a proxima interacao
    dfixi_1 = dfixi;
    
    ## Incrementa o contador
    i = i + 1;      
    
  }
  
  ## Calcula o ultimo polinomio
  ystar = ai + bi*xstar + ci*xstar^2 + di*xstar^3;
  
  ## Retorna o polinomio calculado
  return(ystar);
  
}


PricePU_BR <- function(yield,days,faceValue=100,coupon=0.1,daysYear=252) {
  ## Description: Compute bond price considering work days (year=252)
  ## Arguments:
  ##  yield: YTM
  ##  days: vector days to coupon payment
  ##  faceValue
  ##  coupon rate semmiannual
  ##  wd: logic variable TRUE=consider working days instead of calendar days (BARZIL only)
  ## Return
  ##  price
  price <- 0
  n <- length(days)-1
  i <- 1
  while (i <= n) {
    price <- price + round(round(faceValue*((1+coupon)^0.5-1),5)/(1+yield)^(trunc(days[i]/daysYear*1e14)/1e14),9)
    i <- i+1
  } # end-while
  price <- trunc((price+round(round(faceValue*((1+coupon)^0.5),5)/(1+yield)^(trunc(days[length(days)]/daysYear*1e14)/1e14),9))*1e6)/1e6  
  return(price)
} # end-function


