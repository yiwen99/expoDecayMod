\name{exp_decay_mod}
\alias{exp_decay_mod}
\title{Fitting exponential decay model}
\usage{
exp_decay_mod(params_ini, sigma, t, u)
}
\arguments{
\item{params_ini}{initial values for the parameters to be estimated.}

\item{sigma}{a binary vector taking value 0 and 1 to show whether t is an event time or censor time}

\item{t}{a numeric vector of the censored time or event time}

\item{u}{a numeric vector of the transplant time}
}
\value{
the estimated parameters(beta, gamma)
}
\description{
Find the estimated paramaters of the exponential decay model.
}
\examples{
library(survival)
data(heart)
fit1 <- coxph(Surv(start,stop,event)~(age + surgery)*transplant,data=heart, method = "efron")
summary(fit1)
exp_decay_mod(c(-0.29295,0.001),sigma=heart$event, t=heart$stop, u=heart$start)#initial value -0.29295 is obtained from fit1
}
