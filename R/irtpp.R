#' Estimate a test item parameters according to Item Response Theory.
#' @param dataset The matrix with the responses from the individuals
#' @param model The model used to calibrate the parameters
#' @param dims The dimensions to use on the estimation, remember to use the initial parameters if you want highquality estimation
#' @param initialvalues The matrix with the initial values for the optimization process
#' @param filename Optional argument specifying a CSV file to read instead of a dataset in memory
#' @param output Optional. Additonal arguments that need to be documented by cristian
#' @param loglikflag Optional. Show the loglikelihood at the end of estimation procedure. Also shows AIC and BIC statistic
#' @return The item parameters in a matrix.
#' @export

irtpp <- function(dataset          = NULL,
                  model,
                  dims             = 1   ,
                  initialvalues    = NULL,
                  filename         = NULL,
                  output           = NULL,
                  loglikflag       = F)
{
  if(dims > 1) {}
  else
  {
    mod = irtpp.model(model,asnumber=T);
    ret = uirtestimate(dataset,mod)

    ret$z = ret$z[1:(mod*ncol(dataset))]
    ##print(ret$z)
    ret$z = matrix(ret$z,ncol=mod,byrow=T)
    ##print("mt")
    if(mod == 2)
    {
      ret$z = cbind(ret$z,rep(0,ncol(dataset)))
    }
    if(mod == 1)
    {
      ret$z = cbind(rep(1,ncol(dataset)),ret$z,rep(0,ncol(dataset)))
    }

    ##print(ret$z)
    ##ret$z = parameter.matrix(ret$z, byrow=T)
    ret$z = parameter.list(ret$z)

    ### now do the loglik.
    if(loglikflag)
    {
      z       = ret$z
      theta   = ret$theta
      weights = ret$weights

      ##first calculate prob matrix
      thsum = NULL
      idx   = 1;
      for (th in theta)
      {
        thsum[[idx]] = probability.3pl(z=z,theta=th)
        idx = idx +1;
      }
      thsum
      thmat   = t(matrix(unlist(thsum),ncol = length(theta)))
      i       = 1
      logliks = 0;
      idx     = 1;
      pfrq    = pattern.freqs(data = dataset)
      pfrqs   = pfrq[,ncol(dataset)+1]
      pat     = pfrq[,1:ncol(dataset)]
      pat     = exp(as.matrix(pat)%*%t(log(thmat))+(1-as.matrix(pat))%*%t(log(1-thmat)))

      pwei    = pat%*%weights
      LL      = -sum(log(rep(pwei,pfrqs)))
      print(paste("Loglikelihood : ",LL));
      AIC     = -2*(-LL)+ 2* ncol(dataset)*3;
      BIC     = -2*(-LL)+ log(nrow(dataset)) * ncol(dataset)*3;
      print(paste("AIC : ",AIC));
      print(paste("BIC : ",BIC));
    }
  }
  ret
}

#' Estimate the latent traits of the individuals in a test with some given item parameters
#' @param model The model used to calibrate the parameters
#' @param itempars The item parameters for the model.
#' @param method The method to estimate traits
#' @param dataset The matrix with the responses from the individuals
#' @param probability_matrix The probability matrix in case it does not need to be recalculated
#' @return A list with the patterns and the estimated latent traits
#' @export
individual.traits<-function(model,
                            itempars,
                            method,
                            dataset             = NULL,
                            probability_matrix  = NULL)
{
  model = irtpp.model(model,asnumber=T)
  cuads = as.matrix(read.table(system.file("extdata","Cuads.csv",package="IRTpp"),sep=",",header=T))

  if(is.null(dataset)) { stop("Please provide a dataset or filename") }
  est.method = ifelse(method == "EAP", eapinterface, mapinterface)

  est = est.method(zita_par     = itempars,
                   dat          = dataset,
                   e_model      = model,
                   matrix_flag  = !is.null(probability_matrix),
                   prob_matrix  = probability_matrix
                   )

  est   = individual.traits.aux(dataset=dataset, est=est)
  freqs = pattern.freqs(dataset);
  freqs = freqs[,ncol(freqs)];
  trait = est$trait;
  est   = cbind.data.frame(est$patterns,freqs,trait);
  est
}

individual.traits.aux <- function(dataset, est)
{
  est = list(matrix(est[[1]],ncol=dim(dataset)[[2]],byrow=T),est[[2]])
  names(est) <- c("patterns","trait")
  est
}
