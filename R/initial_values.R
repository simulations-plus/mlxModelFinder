err <-  function(parameter, y, p.ind, N, a) {
  p.ind[,] <- matrix(exp(parameter),nrow=N,ncol=length(parameter),byrow=TRUE)
  f <- as.numeric(lixoftConnectors::computePredictions(p.ind)[[1]])
  if (any(is.nan(f)) | any(is.infinite(f)))
    e <- Inf
  else
    e <- mean((log(f+a)-log(y+a))^2)
  #    e <- mean((f^a-y^a)^2)
  return(e)
}

pk.estim <- function(r, admin) {
  time <- NULL
  g <- lixoftConnectors::getObservationInformation()
  gn <- g$name[[1]]
  gy <- g[[gn]]

  treat <- lixoftConnectors::getTreatmentsInformation()
  if (!is.null(treat$rate))
    treat$tinf <- treat$amount/treat$rate
  else
    treat$tinf <- 0
  abs <- elim <- max1 <- trid <- amount2 <- NULL

  r <- gy
  for (id in r$id) {
    ji <- which(treat$id==id)
    tri <- treat[ji,]
    tri <- tri[order(tri$time),]

    ji <- which(gy$id==id)
    yi <- gy[ji,]
    yi <- yi[yi[['time']]>= min(tri[['time']]),]

    jty1 <- which(tri[['time']]<=min(yi[['time']]))
    tri <- tri[max(jty1):nrow(tri), ]
    jty2 <- which(tri[['time']]>max(yi[['time']]))
    if (length(jty2)>1)
      tri <- tri[1:(min(jty2)-1), ]

    ndi <- nrow(tri)
    tri.inf <- tri[ndi,]
    tri.inf['time'] <- Inf
    tri <- rbind(tri, tri.inf)

    tri1 <- tri[1:2,]
    tri2 <- tri[ndi:(ndi+1),]
    #trid <- rbind(trid, tri1[1,])

    yi1 <- subset(yi, time>=tri1$time[1] & time<tri1$time[2] )
    if (nrow(yi1)>0) {
      yi1[gn] <- yi1[gn]/tri1[['amount']][1]
      yi1['time'] <- yi1['time'] - tri1$time[1]
      if (admin=="oral") {
        j.max1 <- which.max(yi1[[gn]])
        if (j.max1 == 1) {
          yi1 <- rbind(yi1[1,], yi1)
          yi1[1,'time'] <- yi1[1,'time']/4
          yi1[1,gn] <- yi1[1,gn]/2
          j.max1 <- 2
        }
        if (length(j.max1)>0 && j.max1>1) {
          abs <- rbind(abs, yi1[1:(j.max1-1),])
          # if (length(j.max1)>0) {
          #   abs <- rbind(abs, yi1[1:(j.max1),])
          if (j.max1<nrow(yi1))
            max1 <- rbind(max1,yi1[j.max1,])
        }
      } else {
        tinfi <- tri1[1,"tinf"]
        #yi1 <- yi1[yi1$time <= tinfi,]
        yi1$tinf <- tinfi
        abs <- rbind(abs, yi1)
      }
    }

    yi2 <- subset(yi, time>=tri2$time[1] & time<tri2$time[2] )
    if (nrow(yi2)>0) {
      amti2 <- tri2[['amount']][1]
      yi2[gn] <- yi2[gn]/amti2
      amount2 <- c(amount2, amti2)
      yi2['time'] <- yi2['time'] - tri2$time[1]
      yi2['amount'] <- amti2
      j.max2 <- which.max(yi2[[gn]])
      ni <- nrow(yi2)
      if (length(j.max2)>0 && j.max2<ni)
        elim <- rbind(elim, yi2[(j.max2):ni,])
    }
  }
  if (!is.null(abs)) names(abs)[which(names(abs)==gn)] <- "y"
  if (!is.null(elim)) names(elim)[which(names(elim)==gn)] <- "y"
  if (!is.null(max1)) names(max1)[which(names(max1)==gn)] <- "y"
  return(list(abs=abs, elim=elim, max1=max1, amount2=amount2))
}

compute.ini <- function(r, admin, parameter) {

  y <- NULL

  th <- pk.estim(r, admin)
  abs <- th$abs
  elim <- th$elim

  k_ini <- -lm(log(y) ~ time, data=subset(elim, y>0))$coefficients[[2]]

  if (admin=="oral") {
    ymax <- th$max1$y
    tmax <- th$max1$time
    ka_ini <- lm(log(y) ~  time, data=subset(abs, y>0))$coefficients[[2]]
    if (is.na(ka_ini) || ka_ini < 0)
      ka_ini <- 1
    Tk0_ini <- mean(tmax)
    if (ka_ini>0)
      V_ini <- 1/mean(ymax)*ka_ini/abs(ka_ini-k_ini)
    else
      V_ini <- 1/(Tk0_ini*k_ini*mean(ymax))*(1-exp(-k_ini*Tk0_ini))
    Tlag_ini <- Tk0_ini/5
    Mtt_ini <- Tlag_ini
    Ktr_ini <- ka_ini*5
    list.ini <- c(ka=ka_ini, V=V_ini, k=k_ini, Tk0=Tk0_ini, Tlag=Tlag_ini, Mtt=Mtt_ini, Ktr=Ktr_ini)
  } else {
    rV <- (1 - exp(-k_ini*abs$tinf)) /(abs$tinf*k_ini)
    rV[abs$tinf==0] <- 1
    dt <- pmax(abs$time - abs$tinf, 0)
    rV <- rV*exp(-k_ini*dt)
    V_ini <- mean(rV^2)/mean(rV*abs$y)
    list.ini <- c(V=V_ini, k=k_ini)
  }
  Cl_ini <- k_ini*V_ini
  Cmax <- aggregate(elim$y*elim$amount, by=list(elim$id), FUN=max)
  Km_ini <- mean(Cmax[,2])
  Vm_ini <- Cl_ini*(2*Km_ini)
  k12_ini <- k_ini/2
  k21_ini <- k_ini/2
  k13_ini <- k_ini/4
  k31_ini <- k_ini/4

  list.ini <- c(list.ini, Cl=Cl_ini, Km=Km_ini, Vm=Vm_ini)
  list.ini <- c(list.ini, k12=k12_ini, k21=k21_ini, k13=k13_ini, k31=k31_ini)
  list.ini <- c(list.ini, V1=V_ini, Q=k12_ini*V_ini, Q2=k12_ini*V_ini, V2=k12_ini/k21_ini*V_ini,
                Q3=k13_ini*V_ini, V3=k13_ini/k31_ini*V_ini)
  return(list.ini[parameter])
}

pop.opt <- function(p0) {
  #setPKproject(parameter=p0)
  g <- lixoftConnectors::getObservationInformation()
  gn <- g$name[[1]]
  gy <- g[[gn]]
  N <- length(unique(gy[['id']]))
  y <- gy[[gn]]
  if (N>1)
    p.ind <- as.data.frame(t(p0)[rep(1,N),])
  else
    p.ind <- p0
  a <- max(-min(y) + 0.5, 0.5)
  # if ("F" %in% names(p0))  browser()
  # p_ini <<- p.ind
  r <- optim(log(p0), err, y=y, p.ind=p.ind, N=N, a=a)
  return(exp(r$par))
}


pkpopini <- function(project, new.project = NULL, new.dir = NULL) {

  lixoftConnectors::loadProject(project)
  data <- lixoftConnectors::getData()[c('dataFile', 'headerTypes')]
  parameter <- lixoftConnectors::getIndividualParameterModel()$name
  y.name <- lixoftConnectors::getObservationInformation()$name

  if (is.null(data$administration)) {
    if ("ka" %in% parameter | "Tk0" %in% parameter)
      data$administration <- "oral"
    else
      data$administration <- "iv"
  }

  if (!is.null(new.dir) && !dir.exists(new.dir))
    dir.create(new.dir)

  if (length(which(c("Q","Q2") %in% parameter)) > 0)
    param <- "clearance"
  else
    param <- "rate"

  if (is.null(new.project))
    new.project <- paste0("pk_", paste0(parameter, collapse=""), '.mlxtran')
  if (is.null(new.dir))
    new.dir <- "."
  new.project <- file.path(new.dir, new.project)

  g <- lixoftConnectors::getContinuousObservationModel()
  eval(parse(text=paste0('lixoftConnectors::setErrorModel(',names(g$errorModel),'= "combined2")')))

  r <- lixoftConnectors::getInterpretedData()
  pini <- compute.ini(r, data$administration, parameter)

  popt <- pop.opt(pini)
  popt <- popt[parameter]

  pop.ini <- lixoftConnectors::getPopulationParameterInformation()
  j.pop <- which(pop.ini$name %in% paste0(parameter,"_pop"))
  j.pop <- match(paste0(parameter,"_pop"), pop.ini$name[j.pop])
  pop.ini$initialValue[j.pop] <- popt
  lixoftConnectors::setPopulationParameterInformation(pop.ini)

  lixoftConnectors::saveProject(projectFile = new.project)

  return(list(pop.ini=popt, project=new.project, data=data))
}

initialValues <- function(project) {
  tryCatch({
    withCallingHandlers(
      {
        pkpopini(project = project, new.project = "pop_init.mlxtran",
                 new.dir = tempdir())
      },
      message = function(m) {
        invokeRestart("muffleMessage")
        stop(m)
      },
      warning = function(w) {
        invokeRestart("muffleWarning")
        stop(w)
      }
    )
  },
  error = function(e) {
    lixoftConnectors::loadProject(project)
    init_values <- lixoftConnectors::getFixedEffectsByAutoInit()
    lixoftConnectors::setPopulationParameterInformation(init_values)
    lixoftConnectors::saveProject(file.path(tempdir(), "pop_init.mlxtran"))
  })

  return(file.path(tempdir(), "pop_init.mlxtran"))
}

