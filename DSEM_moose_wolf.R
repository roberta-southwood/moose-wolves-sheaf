library(dsem)
library(dynlm)

# dsem gives similar results to a vector autoregressive (VAR) model

# code from Thorson 
vignette("features", package = "dsem")

data(isle_royale)
data = ts( log(isle_royale[,2:3]), start=1959)

#define model 
sem = "
  # Link, lag, param_name
  wolves -> wolves, 1, arW
  moose -> wolves, 1, MtoW
  wolves -> moose, 1, WtoM
  moose -> moose, 1, arM
"
#fit the model 

# initial first without delta0 (to improve starting values)
fit0 = dsem( sem = sem,
             tsdata = data,
             estimate_delta0 = FALSE,
             control = dsem_control(
               quiet = FALSE,
               getsd = FALSE) )
#>   Coefficient_name Number_of_coefficients   Type
#> 1           beta_z                      6  Fixed
#> 2             mu_j                      2 Random

#
parameters = fit0$obj$env$parList()
parameters$delta0_j = rep( 0, ncol(data) )

# Refit with delta0
fit = dsem( sem = sem,
            tsdata = data,
            estimate_delta0 = TRUE,
            control = dsem_control( quiet=TRUE,
                                    parameters = parameters ) )
head(fit)

# dynlm
fm_wolf = dynlm( wolves ~ 1 + L(wolves) + L(moose), data=data )   #
fm_moose = dynlm( moose ~ 1 + L(wolves) + L(moose), data=data )   #

# MARSS
library(MARSS)
z.royale.dat <- t(scale(data.frame(data),center=TRUE,scale=FALSE))
royale.model.1 <- list(
  Z = "identity",
  B = "unconstrained",
  Q = "diagonal and unequal",
  R = "zero",
  U = "zero"
)
kem.1 <- MARSS(z.royale.dat, model = royale.model.1)
#> Success! abstol and log-log tests passed at 16 iterations.
#> Alert: conv.test.slope.tol is 0.5.
#> Test with smaller values (<0.1) to ensure convergence.
#> 
#> MARSS fit is
#> Estimation method: kem 
#> Convergence test: conv.test.slope.tol = 0.5, abstol = 0.001
#> Estimation converged in 16 iterations. 
#> Log-likelihood: -3.21765 
#> AIC: 22.4353   AICc: 23.70964   
#>  
#>                       Estimate
#> B.(1,1)                 0.8669
#> B.(2,1)                -0.1240
#> B.(1,2)                 0.1417
#> B.(2,2)                 0.8176
#> Q.(X.wolves,X.wolves)   0.1341
#> Q.(X.moose,X.moose)     0.0284
#> x0.X.wolves             0.2324
#> x0.X.moose             -0.6476
#> Initial states (x0) defined at t=0
#> 
#> Standard errors have not been calculated. 
#> Use MARSSparamCIs to compute CIs and bias estimates.
SE <- MARSSparamCIs( kem.1 )

# Using VAR package
library(vars)
var = VAR( data, type="const" )

### Compile
m1 = rbind( summary(fm_wolf)$coef[-1,], summary(fm_moose)$coef[-1,] )[,1:2]
m2 = summary(fit$sdrep)[1:4,]
#m2 = cbind( "Estimate"=fit$opt$par, "Std. Error"=fit$sdrep$par.fixed )[1:4,]
m3 = cbind( SE$parMean[c(1,3,2,4)], SE$par.se$B[c(1,3,2,4)] )
colnames(m3) = colnames(m2)
m4 = rbind( summary(var$varresult$wolves)$coef[-3,], summary(var$varresult$moose)$coef[-3,] )[,1:2]

# Bundle
m = rbind(
  data.frame("var"=rownames(m1), m1, "method"="dynlm", "eq"=rep(c("Wolf","Moose"),each=2)),
  data.frame("var"=rownames(m1), m2, "method"="dsem", "eq"=rep(c("Wolf","Moose"),each=2)),
  data.frame("var"=rownames(m1), m3, "method"="MARSS", "eq"=rep(c("Wolf","Moose"),each=2)),
  data.frame("var"=rownames(m1), m4, "method"="vars", "eq"=rep(c("Wolf","Moose"),each=2))
)
#knitr::kable( m1, digits=3)
#knitr::kable( m2, digits=3)

m = cbind(m, "lower"=m$Estimate-m$Std..Error, "upper"=m$Estimate+m$Std..Error )

# ggplot estimates ... interaction(x,y) causes an error sometimes
library(ggplot2)
library(ggpubr)
library(ggraph)
longform = reshape( isle_royale, idvar = "year", direction="long", varying=list(2:3), v.names="abundance", timevar="species", times=c("wolves","moose") )
p1 = ggplot( data=longform, aes(x=year, y=abundance) ) +
  facet_grid( rows=vars(species), scales="free" ) +
  geom_point( )

p2 = ggplot(data=m, aes(x=interaction(var,eq), y=Estimate, color=method)) +
  geom_point( position=position_dodge(0.9) ) +
  geom_errorbar( aes(ymax=as.numeric(upper),ymin=as.numeric(lower)),
                 width=0.25, position=position_dodge(0.9))  #
p3 = plot( as_fitted_DAG(fit, lag=1), rotation=0 ) +
  geom_edge_loop( aes( label=round(weight,2), direction=0)) + #arrow=arrow(), , angle_calc="along", label_dodge=grid::unit(10,"points") )
  expand_limits(x = c(-0.1,0) )

ggarrange( p1, p2, p3,
           labels = c("Time-series data", "Estimated effects", "Fitted path digram"),
           ncol = 1, nrow = 3)

# Calculate total effects
effect = total_effect( fit )

# Plot total effect
ggplot( effect) + 
  geom_bar( aes(lag, total_effect, fill=lag), stat='identity', col='black', position='dodge' ) +
  facet_grid( from ~ to  )




data("isleRoyal")

log_isle_royale <- log(isle_royale)

write.csv(log_isle_royale, "Downloads/research/log_isle_royal_data.csv")
write.csv()


log_dsem_values <- (fit$internal$parhat$x_tj)

write.csv(log_dsem_values, "Downloads/research/wolf_moose_dsem.csv")

head(fit$internal$parhat$x_tj)







