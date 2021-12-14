\name{exp_decay_mod}
\alias{exp_decay_mod}
\title{Fitting exponential decay model}
\usage{
exp_decay_mod(params_ini, sigma, t, u)
}
\arguments{
\item{params_ini}{initial values for parameters /beta and /gamma.}

\item{sigma}{a binary vector taking value 0 when t is a censoring time and 1 when t is an event time}

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
# For large dataset, the analysis might not be able to run on local computers. Here is our example analysis supported by Google Collab Notebook: https://colab.research.google.com/drive/1r99N7ws5RECkjBDUt6chTz2cebkdvrtd#scrollTo=mEj1bfQpbSLI
}
