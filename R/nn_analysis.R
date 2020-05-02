library(hadron)
library(Raff)
library(pracma)
library(binhf)
library(boot)
library(DescTools)
library(stringr)
source('/qbigwork2/beilschmidt/code/R/get_sl.R')
source('/qbigwork2/beilschmidt/code/R/R_gamma.R')
#options(warn = -1)#, error=recover)
trace <- function(A) {
  n <- dim(A)[1] # get dimension of matrix
  tr <- 0 # initialize trace value
  
  # Loop over the diagonal elements of the supplied matrix and add the element to tr
  for (k in 1:n) {
    l <- A[k,k]
    tr <- tr + l
  }
  return(tr[[1]])
}

projector <- function(sgn, gamma0){
  p <- array(0, dim = c(4,4))
  for (a in c(1:4)) {
    for (b in c(1:4)) {
      if (a == b) {
        p[a,b] <- 0.5
      } else {
        p[a,b] <- 0.5 * sgn * gamma0[a,b]
      }
    }
  }
  return(p)

}



read_dia <- function(key, nrDiagram=c(1:4), file_str, src_str_alt, .gi= "Gi_Cg5", .gf = "Gf_Cg5", .mom_tag){

    diagrams <- list()
    corrKey <- key
    for (t in nrDiagram) {
          t_name <-  paste("t",toString(t), sep = "")
          key_str <- paste("/",corrKey,"/", src_str_alt,"/",.gf,"/Gc_id/",.gi,"/QX0_QY0_QZ0/",t_name,"/", .mom_tag , sep = "")
          d <- aff_read_key(file_str, key_str, (4*4*time))
          #diagrams[t_name] <- raw_cf_data(cf = raw_cf_meta(Time=T, dim=c(4,4)),
          #                                  data = aperm(array(d, dim = c(4, 4, 64)), c(3,1,2)))
          diagrams[[t_name]] <- array(d, dim = c(4, 4, time))
     }
    return(Reduce('+', diagrams))    
}
  pathToData<-function(letter){
      return(sprintf("/hiskp4/petschlies/nucleon-ff/cA211%s.30.32/J125-k0p2/", letter))
  }

read_dia_list <- function(key, nrDiagram=c(1:4), file_str, src_str_alt, .gi= "Gi_Cg5", .gf = "Gf_Cg5", .mom_tag, two_pt=FALSE){

    key_str <- c()
    corrKey <- key
    pre <- "t"
    if(two_pt) pre <- "n"
    for (t in nrDiagram) {
          t_name <-  paste(pre,toString(t), sep = "")
          if(two_pt){ str <- paste("/","N-N","/", src_str_alt,"/",.gi,"/",.gf,"/",t_name,"/", .mom_tag, sep = "")}
          else { str <- paste("/",corrKey,"/", src_str_alt,"/",.gf,"/Gc_id/",.gi,"/QX0_QY0_QZ0/",t_name,"/", .mom_tag , sep = "")}
          key_str <-c(key_str, str )
          #d <- aff_read_key(file_str, key_str, (4*4*time))
          #diagrams[t_name] <- raw_cf_data(cf = raw_cf_meta(Time=T, dim=c(4,4)),
          #                                  data = aperm(array(d, dim = c(4, 4, 64)), c(3,1,2)))
          #diagrams[[t_name]] <- array(d, dim = c(4, 4, time))
     }
    d <- aff_read_key_list(file_str, key_str, (4*4*time))
    d <- array(d, dim=c(4, 4, time, length(nrDiagram)))
    return(apply(d, MARGIN=c(1,2,3), sum))    
}

  pathToData<-function(letter){
      return(sprintf("/hiskp4/petschlies/nucleon-ff/cA211%s.30.32/J125-k0p2/", letter))
  }

getSrcList <- function(letter){

  if(strcmp(letter, "a")){

        end_str = ".all.tab"
    } else {
        end_str = ".tab"
    }
  return( source_list <- get_source_coords_table(f = paste(pathToData(letter),"source_coords.cA211",letter,".30.32.nsrc16",end_str, sep = "")))
}
###############################################################
#main function for analyzing correlator function 
###############################################################

#' analyse Function
#'
#' This function allows you to read a specific correlator from an aff-file.
#' @param 
#' @keywords
#' @export
#' @examples
#' cf <- analyze()

analyse <- function(corrKey = "N-N",T=time, Lx=time/2, Ly=time/2, Lz=time/2, n_src = 16,n_conf = 1224, path_letter = "b", gi= "Gi_Cg5", gf = "Gf_Cg5", step = 24, mom_tag = "px00py00pz00", conf_start=0, meanOverSrc=TRUE){

    
    num_conf <- (n_conf-conf_start)/step+1
    cf <- array(0, dim = c(T, (nchar(path_letter)*num_conf)))
    l <- 0

    letter_vec <- strsplit(path_letter,"")[[1]]
    for(letter in letter_vec){ 

    l <- l+1
    cor <- array(0, dim = c(T, n_src, (num_conf)))
        cor_m <- array(0, dim = c(T, n_src, (num_conf)))
        cor_mean_src <- array(0, dim = c(T, (num_conf)))
        #load gamma matrices
        g <- set_gamma_basis_all(g=gamma_basis$tmlqcd)
        #load Cxgamma, but without i factor and multiply with i
        Cg <- set_Cgamma_basis_matching( g=gamma_basis$tmlqcd  )
  
        # generate projektor matrix
        proj_p <- projector(1,g$`0`)
        proj_m <- projector(-1,g$`0`)
        # Source coord list
        #"/hiskp4/petschlies/nucleon-ff/cA211a.30.32"
        source_list <- getSrcList(letter) 
        print(sprintf("%s, T=%i, n_src=%i, n_conf=%i, Gi = %s, Gf = %s, mom = %s", corrKey, T, n_src, (num_conf),gi, gf, mom_tag))
        for (conf in seq(conf_start, n_conf, step)){
            for (src in c(1:n_src)) {
      
                src_str <- get_source_coords_tag(c = conf, i=src, sl = source_list)
     
                # Open File
                file_str <- paste(pathToData(letter),"njn_fht.",formatC(conf, width = 4, flag = "0"),".", src_str,".aff", sep="")
                src_str_alt <- toupper(src_str)
              
                for (i in c(2:4)) {
                    tmp_str <- unlist(strsplit(src_str_alt, ""))
                    ind <- grep("[A-Z]", tmp_str)[i]
                    src_str_alt <- paste0(c(tmp_str[1:(ind-1)], "_", tmp_str[ind:(nchar(src_str_alt))]), collapse="")
                }
              
                if(strcmp(corrKey, "P-P")){
                    #key_str <- paste("/","N-N","/", src_str_alt,"/",gi,"/",gf,"/n1/", mom_tag, sep = "")
                    #d <- aff_read_key(file_str, key_str, (4*4*T))
                    # create SpinxSpinxTime matrix n1
                    #n1 <- array(d, dim = c(4, 4, T))
                    #key_str <- paste("/","N-N","/", src_str_alt,"/",gi,"/",gf,"/n2/", mom_tag, sep = "")
                    #d <- aff_read_key(file_str, key_str, (4*4*T))

                    # create SpinxSpinxTime matrix n2
                    #n2 <- array(d, dim = c(4, 4, T))

                    nn <-read_dia_list("N-N", nrDiagram=c(1:2), file_str, src_str_alt, .mom_tag=mom_tag, two_pt=TRUE) #n1 + n2
                
                } else if (strcmp(corrKey, "N-N") & strcmp(letter, "b")) {
                    
                    #key_str <- paste("/","N-N","/", src_str_alt,"/",gi,"/",gf,"/n3/", mom_tag, sep = "")
                    #d <- aff_read_key(file_str, key_str, (4*4*T))
                    ## create SpinxSpinxTime matrix n1
                    #n1 <- array(d, dim = c(4, 4, T))
                    #key_str <- paste("/","N-N","/", src_str_alt,"/",gi,"/",gf,"/n4/", mom_tag, sep = "")
                    #d <- aff_read_key(file_str, key_str, (4*4*T))

                    ## create SpinxSpinxTime matrix n2
                    #n2 <- array(d, dim = c(4, 4, T))


                    #nn <- n1 + n2
                    nn <-read_dia_list("N-N", nrDiagram=c(3:4), file_str, src_str_alt, .mom_tag=mom_tag, two_pt=TRUE)
                } else if (strcmp(corrKey, "P-ubGu-P")) {
                    
                    key <- "N-ubGu-N"
                    nn <- read_dia_list(key, nrDiagram =c(1:4) , file_str, src_str_alt, .mom_tag=mom_tag, .gi=gi, .gf=gf)
                } else if (strcmp(corrKey, "P-dbGd-P")) {
                    key <- "N-dbGd-N"
                    nn <- read_dia_list(key, nrDiagram =c(1:2) , file_str, src_str_alt, .mom_tag=mom_tag, .gi=gi, .gf=gf)

                } else if (strcmp(corrKey, "P-J-P") & strcmp(letter, "a")) {
                    nn <- read_dia_list("N-ubGu-N", nrDiagram =c(1:4) , file_str, src_str_alt, .mom_tag=mom_tag, .gi=gi, .gf=gf)
                    nn <- nn + read_dia_list("N-dbGd-N", nrDiagram =c(1:2) , file_str, src_str_alt, .mom_tag=mom_tag, .gi=gi, .gf=gf)

                } else if (strcmp(corrKey, "P-J-P") & strcmp(letter, "b")){
                    
                    nn <- read_dia_list("N-qbGq-N", nrDiagram =c(1:6) , file_str, src_str_alt, .mom_tag=mom_tag, .gi=gi, .gf=gf)

                }  else if (strcmp(corrKey, "N-J-N") & strcmp(letter, "b")){
                    
                    nn <- read_dia_list("N-qbGq-N", nrDiagram =c(7:12) , file_str, src_str_alt, .mom_tag=mom_tag, .gi=gi, .gf=gf)
                } else {
                    
                    stop("No valid corrKey and letter combination given to analyse!")
                }


                #is.raw_cf(nn)
                #get_plotdata_raw_cf(nn, "both", TRUE, TRUE) 
                # get soruce time
                src_mat <- matrix(unlist(strsplit(src_str_alt, "_")))
                t_src <- strtoi(substring(src_mat[1], 2, nchar(src_mat[1])))
                x_src <- strtoi(substring(src_mat[2], 2, nchar(src_mat[2])))
                y_src <- strtoi(substring(src_mat[3], 2, nchar(src_mat[3])))
                z_src <- strtoi(substring(src_mat[4], 2, nchar(src_mat[4])))
                
                src_vec <- c( x_src, y_src, z_src)#t_src,

                mom_mat <- matrix(unlist(strsplit(mom_tag, "p")))
                x_mom <- strtoi(substring(mom_mat[2], 2, nchar(mom_mat[2])))
                y_mom <- strtoi(substring(mom_mat[3], 2, nchar(mom_mat[3])))
                z_mom <- strtoi(substring(mom_mat[4], 2, nchar(mom_mat[4])))
                
                mom_vec <-2*pi*c(x_mom/Lx, y_mom/Ly, z_mom/Lz)
                # add phasefactor
                shifted <- array(0, dim=c(4,4,T))
                for (t in c(1:T)) {
                    del_t <- (T + t-1 - t_src) %% T
                    phasefactor <- exp(1i * 3 * pi * del_t / T)
                    nn[,,t] <- nn[,,t]*phasefactor
                    shifted[,,(del_t+1)] <- nn[,,t]
                    # multiply projection matrix and take trace
                    #cor[t,src,((conf-conf_start)/step+1)] <- sum(diag(proj_p %*% nn[,,t]))
                    #cor_m[t,src,((conf-conf_start)/step+1)] <- sum(diag(proj_m %*% nn[,,t]))
                } # end time loop
                
                    impulsfactor <- exp(-1i * dot(mom_vec ,src_vec))
                    nn <- shifted*impulsfactor
                cor[,src,((conf-conf_start)/step+1)] <- apply(nn, MARGIN=3, FUN=function(mat, proj){sum(diag(proj %*% mat)) }, proj=proj_p)
                cor_m[,src,((conf-conf_start)/step+1)] <- apply(nn, MARGIN=3, FUN=function(mat, proj){sum(diag(proj %*% mat)) }, proj=proj_m)

                cor_m[2:T,src,((conf-conf_start)/step+1)] <- Rev(cor_m[2:T,src,((conf-conf_start)/step+1)])
              
                # shift source location to 0
                #cor[,src,((conf-conf_start)/step+1)] <- shift(cor[,src,((conf-conf_start)/step+1)], t_src, dir = "left")
                #cor_m[,src,((conf-conf_start)/step+1)] <- shift(cor_m[,src,((conf-conf_start)/step+1)], t_src, dir = "left")
        } # end source loop
    
   
  
    } # end loop conf
 
    #
    #  print("Not symmetrised positive parity:")
    #  print(apply(Re(apply(cor, MARGIN=c(1,3), mean)), MARGIN=1, mean))
    #
    #  print("Not symmetrised negative parity:")
    #  print(apply(Re(apply(cor_m, MARGIN=c(1,3), mean)), MARGIN=1, mean))

    cor <- (cor - cor_m)/2
    if(!meanOverSrc){
        return(list("cor"=cor, "cor_m"=cor_m))
    }
    cor_mean_src <- apply(cor, MARGIN=c(1,3), mean)
    cf[,((l-1)*num_conf+1):(l*num_conf)] <- cor_mean_src
 } #end loop path_letter
  
  cf <- cf_meta(.cf=cf_orig(cf=Re(t(cf[1:(T/2+1),]))), T=T, symmetrised=TRUE)
  cf <- bootstrap.cf(cf)

  return(cf)
}

fit.constant <- function(M, y) {
        res <- list()
    if(is.matrix(M)){ # this is the covariance case
                m.eff <- sum(M %*% y)/sum(M)
            res$value <- (y-m.eff) %*% M %*% (y-m.eff)
                }else{ # this is the uncorrelated case
                            m.eff <- sum(M*y)/sum(M)
                    res$value <- sum(M*(y-m.eff)^2)
                        }
        res$par <- c(m.eff)
            return(res)
}


ratio <- function( c3pt, indices, c2pt, ttau=1, T ){
  c3 <- apply ( c3pt[indices,], c(2), mean )
  c2 <- apply(c2pt[indices,], c(2), mean)
  idt1 <- ( 0 : ( T - ttau ) ) + 1
  idt2 <- idt1 + ttau

  return ( ( Re(c3[idt1]) / Re(c2[idt1]) - Re(c3[idt2])/Re(c2[idt2]) ) / ttau )
}



bootstrap.ratio <- function(c2, c3, R=5000, tau =1){
    if(length(c2$cf) != length(c3$cf)){
    ind <- c((length(c3$cf[,1])-length(c2$cf[,1])):-1) 
    ind2 <- c(1:length(c3$cf[,1]))
    } else {
        ind <- c((length(c3$cf[,1])/2+1):length(c3$cf[,1]))
        ind2 <- ind
    }
    r <- boot(data = c3$cf[ind2,], R=R, statistic = ratio, c2pt=c2$cf[ind,], ttau=tau, T=c2$T)
    #write.table( df, file="m_boot.tmp", col.names=F, row.names=F )
    return( invisible ( r ) )
}

#' summary_ratiofit Function 
#'
#' This function allows you to summarize the ratio fit.
#' @param 
#' @keywords
#' @export
#' @examples
#' 
summary_ratiofit <- function(object, ..., verbose = FALSE) {
      ratio <- object
    cat("\n ** Result of ratio analysis **\n\n")

    #cat("no. measurements\t=\t", ratio$N, "\n")
    #cat("type\t=\t", ratio$type, "\n")
    cat("boot.R\t=\t", ratio$boot.R, "\n")
    cat("boot.l\t=\t", ratio$boot.l, "\n")
    cat("Time extend\t=\t", ratio$Time, "\n")
    cat("NA count in fitted bootstrap samples:\t", length(which(is.na(ratio$cf[,ratio$ii]))),
        "(",100*length(which(is.na(ratio$cf[,ratio$ii])))/ length(ratio$cf[,ratio$ii]), "%)\n")
    #cat("NAs replaced in fit:", ratio$ratiofit$replace.na, "\n")
    cat("time range from", ratio$ratiofit$t1, " to ", ratio$ratiofit$t2, "\n")
    cat("tau =", ratio$tau, "\n")
    #cat("No correlation functions", ratio$nrObs, "\n")
    #if(verbose) {
    #    cat("values with errors:\n\n")
    #    print(data.frame(t= ratio$t, m = ratio$t0, dm = ratio$se))
    #    }
    #cat("correlated fit\t=\t", ratio$ratiofit$useCov, "\n")
    #cat("m\t=\t", ratio$ratiofit$t0[1], "\n")
    #cat("dm\t=\t", ratio$ratiofit$se[1], "\n")
    cat("R\t=\t", ratio$ratiofit$t0[1], "\n")
    cat("dR\t=\t", ratio$ratiofit$se[1], "\n")
    cat("chisqr\t=\t", ratio$ratiofit$chisqr, "\n")
    cat("dof\t=\t", ratio$ratiofit$dof, "\n")
    cat("chisqr/dof=\t",
                        ratio$ratiofit$chisqr/ratio$ratiofit$dof, "\n")
    cat("Quality of the fit (p-value):",   ratio$ratiofit$Qval, "\n")

}

#' plot_ratio Function
#'
#' This function allows you to plot the ratio of two correlator functions.
#' @param 
#' @keywords
#' @export
#' @examples
#' 

plot_ratio <- function(c2, c3 ,t1, t2, add=FALSE, .col="black", .tau=1,  ...){

    df <- bootstrap.ratio(c2, c3, R=5000, tau=.tau)
    std = apply(df$t, c(2), sd )
    std = std/sqrt(length(df$t[,1]))
    
    if(!add){
   plotwitherror(x=c(1:(c2$T)) , y=Re(df$t0[1:(c2$T)]), dy=Re(std[1:(c2$T)]) , main =paste0( "N-J-N linear response of effective mass to external bilinear current, tau=", .tau),ylab = "g_00",xlab = "t/a", xlim=c(0,30), col=.col,  pch=1,cex=0.8, cex.main=1,lwd = 0.3, frame.plot=FALSE, ...) 
    } else {

   pointswitherror(x=c(1:(c2$T)) , y=Re(df$t0[1:(c2$T)]), dy=Re(std[1:(c2$T)]) , main =paste0( "N-J-N linear response of effective mass to external bilinear current, tau=", .tau),ylab = mtext(dM[eff]/ dlambda),xlab = "t/a", xlim=c(0,30), col=.col,  pch=1,cex=0.8, cex.main=1,lwd = 0.3, frame.plot=FALSE, ...) 
    }
  ratio <- cf_meta(.cf=cf_orig(cf=df$t), T=c2$T, symmetrised=TRUE)
  ratio <- bootstrap.cf(ratio)
  ratio$se <- std
  M <- try(invertCovMatrix(ratio$cf[,t1:t2], boot.samples=TRUE), silent=TRUE)
  if(inherits(M, "try-error")) {
    warning("Inversion of variance covariance matrix failed, continuing with uncorrelated chi^2\n")
    M <- diag(1/ratio$tsboot.se[t1:t2]^2)
  }
  ratio$ii <- c(t1:t2)
  ratio$dof <- length(ratio$ii)-1
  ratiofit.tsboot <- array(0, dim=c(ratio$boot.R,2))
  for(i in c(1:ratio$boot.R)) {
            opt <- fit.constant(M=M, y = ratio$cf[i,t1:t2])
            ratiofit.tsboot[i, 1] <- opt$par[1]
            ratiofit.tsboot[i, 2] <- opt$value
            }
  ratio$ratiofit.tsboot <- ratiofit.tsboot
  ratio$ratiofit <- fit.constant(M=M, y = ratio$cf0[t1:t2])
  ratio$ratiofit$t1 <- t1
  ratio$ratiofit$t2 <- t2
  ratio$tau <- .tau
  ratio$ratiofit$t <- ratio$ratiofit.tsboot 
  ratio$ratiofit$t0 <- c(ratio$ratiofit$par[1], ratio$ratiofit$value)
  ratio$ratiofit$se <- sd(ratiofit.tsboot[c(1:(dim(ratiofit.tsboot)[1]-1)),1] )
  ratio$ratiofit$cf <- ratio$cf
  ratio$ratiofit$ii <- ratio$ii
  ratio$ratiofit$dof <- ratio$dof
  ratio$chisqr <- ratio$ratiofit$value
  ratio$ratiofit$chisqr <- ratio$ratiofit$value
  ratio$Qval <- 1-pchisq(ratio$chisqr, ratio$dof)
  lines(x=c(t1,t2),
                   y=c(ratio$ratiofit$par[1],ratio$ratiofit$par[1]),
                             col=.col,
                             lwd=1.3)
  summary_ratiofit(ratio)
        pcol <- col2rgb(.col,alpha=TRUE)/255                                                                                                   
        pcol[4] <- 0.65
              pcol <- rgb(red=pcol[1],green=pcol[2],blue=pcol[3],alpha=pcol[4])
              rect(xleft=ratio$ratiofit$t1, ybottom=ratio$ratiofit$t0[1]-ratio$ratiofit$se[1],
                              xright=ratio$ratiofit$t2, ytop=ratio$ratiofit$t0[1]+ratio$ratiofit$se[1],
                                         col=pcol,
                                         border=NA)
  return(ratio)

#if(doTau){
#   data <- list()
#   data$tau1 <- df
#    for (ttau in c(2:par.tau)) { # 2:tau need change
#                   
#        df <- bootstrap.ratio(c2,c3, R = 5000, tau = ttau)
#        points(x=df$t[1:(c2$T-1)]+1, y=Re(df$m[1:(c2$T-1)]), pch=ttau,cex=0.8, col=colVec[ttau], lwd = 0.3)
#        arrows((df$t[1:(c2$T-1)]+1), Re(df$m[1:(c2$T-1)]-df$std[1:(c2$T-1)]),(df$t[1:(c2$T-1)]+1), Re(df$m[1:(c2$T-1)]+df$std[1:(c2$T-1)]), length=0.05, angle=90, code=3, lwd = .3)
#   data[[sprintf("tau%i", ttau)]] <- df
#    }
#	legend("topleft",c("tau = 1","tau = 2","tau = 3","tau = 4","tau = 5","tau = 6","tau = 7","tau = 8"),cex=.8,col=colVec,pch=c(1:8))
#    return(invisible(data))
#}
}



pointswitherror <- function(x, y, dy, col="black", ...) {
      points(x=x, y=y, col=col, ...)
      arrows(x0=x, y0=y-dy, x1=x, y1=y+dy, length=0.05,
                  angle=90, code=3, col=col, lw = 0.3)
}


plot_tau <- function(data, n_points = 32, n_tau = 8){
    tau_array <- array(0, dim=c(n_points,n_tau, 2))
    for(tau in c(1:n_tau)){ #change 2 back to 1
         temp <- as.matrix(data[[sprintf("tau%i", tau)]]["m"])
         tau_array[,tau, 1] <- temp[1:n_points]
         temp <- as.matrix(data[[sprintf("tau%i", tau)]]["std"])
         tau_array[,tau, 2] <- temp[1:n_points]

    }
    if(gi== "Gi_Cg5"){
        ylim <- c(-40, 30)
    }else{
        ylim <- c(-30, 30)
    }
    plotwitherror(x=c(1:n_tau),y= tau_array[1,,1], ylim=ylim, dy = tau_array[1,,2], ylab = "FHT-ratio",  xlab = "tau")
    for(time in c(2:n_points)){
        pointswitherror( x=c(1:n_tau),y= tau_array[time,,1], dy = tau_array[time,,2])
    }
}

bootmean <- function(vec, i) mean(vec[i])

plot_ratiofit_vs_tau <- function(part="P", .gi = "Gi_Cg5", .gf = "Gf_Cg5", mom="0", .col="BLUE"){
    
    p_tag <- paste0("p_tot", mom)

    data <- list()
    tau_vec <- c(1:par.tau)
    data$tau <- tau_vec
    for(tau in tau_vec){
        
        tau_str <- paste0("tau", tau)
        ratio <- ratios[[part]][[p_tag]][[.gi]][[.gf]][["ratiofit"]][[paste0("t1_",par.t1)]][[paste0("t2_",par.t2)]][[tau_str]]
        data$val <- c( data$val, ratio$ratiofit$par)
        data$se <- c(data$se , ratio$ratiofit$se)

        if(exists("mat")){
            mat <- cbind(mat, ratio$cf0[1:20])
        } else {
            mat <- matrix(ratio$cf0[1:20])
        }

    }
    
    pdf(sprintf("output/NJN_%s_%s_%s_%s_%i_%i_to_%i_ratio_vs_tau.pdf",part, .gi, .gf, p_tag, test, par.t1, par.t2))
    plotwitherror(x=data$tau, y=data$val, ylim=c(-20, -5), dy = data$se, main = paste0("ratiofit vs tau, particle=", part, ", ", .gi, ", ", .gf, ", ",p_tag, ",\nt1=", par.t1, ", t2=", par.t2 ), ylab="Ratiofit Data",  xlab = "tau")

    
    #Cinv <- diag(1/data$se^2)
    C <- cov(mat)
    Cinv <- solve(C)
    Z <- matrix(rep(1, par.tau))
    M <- t(Z) %*% Cinv %*% Z
    Mchol <- chol(M)
    opt <- list()
    opt$sd <- chol2inv(Mchol)
    opt$se <- opt$sd/sqrt(par.tau)
    opt$par <- opt$sd %*% t(Z) %*% Cinv %*% data$val

  lines(x=c(1,par.tau),
                   y=c(opt$par[1],opt$par[1]),
                             col=.col,
                             lwd=1.3)
    #bs <- boot(data$val, bootmean, R=4000, stype="i")
    #opt$se <- sd(bs$t)/sqrt(par.tau)
        pcol <- col2rgb(.col,alpha=TRUE)/255                                                                                                   
        pcol[4] <- 0.65
              pcol <- rgb(red=pcol[1],green=pcol[2],blue=pcol[3],alpha=pcol[4])
              rect(xleft=1, ybottom=opt$par[1]-opt$se[1],
                              xright=par.tau, ytop=opt$par[1]+opt$se[1],
                                         col=pcol,
                                         border=NA)
   
    legend("topleft", c(paste0("fit value: ", opt$par[1], " pm ", opt$se)))
dev.off()

}



calc_cf <- function(.gi, .gf, .mom_tag, test, letter, particle = "P" ){
    str_2pt <- paste(particle, "-", particle, sep="")
    str_3pt <- paste(particle, "-J-", particle, sep="") 
    c2<- analyse(corrKey = str_2pt, gi= .gi, gf = .gf, n_conf=max_conf, step = stepwidth, mom_tag=.mom_tag, path_letter = letter, conf_start = min_conf)
    if(strcmp(.gi, "Gi_Cg5gt") || strcmp(.gf, "Gf_Cg5gt")){
        l <- "b"
    } else {
        l <- letter
    }
    c3 <- analyse(corrKey = str_3pt, gi=.gi, gf = .gf, n_conf=max_conf, step = stepwidth, mom_tag = .mom_tag, path_letter = l, conf_start = min_conf)

    return(invisible(list("c2"=c2, "c3"=c3)))
}

#' calc_all Function to calculate all 2pt and 3pt functions for every momemtum in mom_tag and and source sink matrix combination in gi.
#'
#' Uses calc all 2pt cfs 
#' @param 
#' @keywords
#' @export
#' @examples

calc_all <- function(test, letter, particle="P"){

    if(strcmp(particle, "N") & !strcmp(letter, "b") ){

        stop("Neutron data only in b-stream!")

    }
        for(p_tag in mom_tag){
            num_i <- str_count(p_tag, "1")
            factor <-1
            if (num_i==1){
                factor <- 1./6
            }else if (num_i==2){
                factor<- 1./12.
            } else if (num_i==3){
                factor <- 1./8
            }

            if(is.null(cf_2pt[[particle]])){
                cf_2pt[[particle]] <<- list()
                cf_3pt[[particle]] <<- list()


            }
            p_tot_tag <- sprintf("p_tot%i", num_i)
            if(is.null(cf_2pt[[particle]][[p_tot_tag]])){
                cf_2pt[[particle]][[p_tot_tag]] <<- list()
                cf_3pt[[particle]][[p_tot_tag]] <<- list()


            }
            for( g in gi){
           # d[[sprintf("%s_%s_%s", g[1], g[2], p_tag)]] <- calc_cf(g[1], g[2], p_tag, test)
            cf <- calc_cf(g[1], g[2], p_tag, test, letter, particle)
            if(is.null(cf_2pt[[particle]][[p_tot_tag]][[g[1]]])){
                
                cf_2pt[[particle]][[p_tot_tag]][[g[1]]] <<- list()
                cf_3pt[[particle]][[p_tot_tag]][[g[1]]] <<- list()

            }
            if(is.null(cf_2pt[[particle]][[p_tot_tag]][[g[1]]][[g[2]]])){
                cf_2pt[[particle]][[p_tot_tag]][[g[1]]][[g[2]]] <<- mul.cf(cf$c2 , factor)
                cf_3pt[[particle]][[p_tot_tag]][[g[1]]][[g[2]]] <<- mul.cf(cf$c3 , factor)

            } else {
                cf_2pt[[particle]][[p_tot_tag]][[g[1]]][[g[2]]] <<- add.cf(cf_2pt[[particle]][[p_tot_tag]][[g[1]]][[g[2]]], cf$c2, b=factor)
                cf_3pt[[particle]][[p_tot_tag]][[g[1]]][[g[2]]] <<- add.cf(cf_3pt[[particle]][[p_tot_tag]][[g[1]]][[g[2]]], cf$c3, b=factor)

            }
        }
    }
}

#' gevp Function for specific 2pt cf
#'
#' Uses GEVP on 2pt cfs 
#' @param 
#' @keywords
#' @export
#' @examples

gevp_2pt <- function(p_tag = "p_tot0", g1 = "Cg5", g2= "Cg5gt", particle="P" ){

        g1i <- paste("Gi_",g1, sep="")
        g1f <- paste("Gf_", g1, sep="")
        g2i <- paste("Gi_",g2, sep="")
        g2f <- paste("Gf_", g2, sep="")
       mean_offdia <- add.cf(cf_2pt[[particle]][[p_tag]][[g1i]][[g2f]], cf_2pt[[particle]][[p_tag]][[g2i]][[g1f]], a=0.5, b=0.5)
        gevp_cf <- c(cf_2pt[[particle]][[p_tag]][[g1i]][[g1f]], mean_offdia)
        gevp_cf <- c(gevp_cf, mean_offdia)
        gevp_cf <- c(gevp_cf, cf_2pt[[particle]][[p_tag]][[g2i]][[g2f]])
        
        gevp_cf <- bootstrap.cf(gevp_cf)
        gevp_cf <- bootstrap.gevp(gevp_cf, element.order=c(1:(gevp_cf$nrObs)))
        
        amp <- gevp2amplitude(gevp_cf, ratios[[particle]][[p_tag]][[g1i]][[g1f]]$massfit)#, type="log")
        plot.gevp.amplitude(amp, main=paste("Gevp Amplitude Plot of ",g1, " and ", g2, " at ", p_tag, " for particle=", particle, sep="" ), ylab="gevp", xlab="t")
        return(invisible(list("amplitude"=amp, "gevp"=gevp_cf)))
}


calc_divideby <- function(factor=1/2){


        for(p_tag in mom_tag){
            num_i <- str_count(p_tag, "1")

            p_tot_tag <- sprintf("p_tot%i", num_i)
            for( g in gi){
                cf_2pt[[particle]][[p_tot_tag]][[g[1]]][[g[2]]] <<- mul.cf(cf_2pt[[particle]][[p_tot_tag]][[g[1]]][[g[2]]] , factor)
                cf_3pt[[particle]][[p_tot_tag]][[g[1]]][[g[2]]] <<- mul.cf(cf_3pt[[particle]][[p_tot_tag]][[g[1]]][[g[2]]] , factor)
            }
        }

}
#' plot Function for 2 correlator
#'
#' This function allows you to plot one specific correlator (2pt and 3pt) combination.
#' @param 
#' @keywords
#' @export
#' @examples

plot_2cf <- function(cf1, cf2, tag1="Cg5", tag2="Cg5gt", col1="blue", col2="red", mom="0", particle="P", cor_type="2pt", ...){
    c1 <- bootstrap.cf(cf1)
    c2 <- bootstrap.cf(cf2)
    t = c(0:(cf1$Time/2))

    title <- paste0(cor_type," Cor. on cA211a/b.30.32, mom=", mom, ",\nparticletype=", particle)
    plotwitherror(x = t, y = c1$cf0, dy = c1$tsboot.se, col=col1, log="y", main=title, xlab="t/a", ylab="C(t)",   ...)
    pointswitherror(x = t, y = c2$cf0, dy = c2$tsboot.se, col=col2, log="y",  ...)
    legend("topright", legend=c(tag1, tag2),
                  col=c(col1, col2), bty="n", pch=c(21,22))
}



#' plot Function for single combination
#'
#' This function allows you to plot one specific correlator (2pt and 3pt) combination.
#' @param 
#' @keywords
#' @export
#' @examples

plot_comb <- function(g.i, g.f, mom.tag, bool, particle="P"){
    plotInfo <- list() 
    if(calcAll){
        num_i <- str_count(mom.tag, "1")
        p_tot_tag <- sprintf("p_tot%i", num_i)
        c2 <- bootstrap.cf(cf_2pt[[particle]][[p_tot_tag]][[g.i]][[g.f]])
        c3 <- bootstrap.cf(cf_3pt[[particle]][[p_tot_tag]][[g.i]][[g.f]])
    } else {
        c2 <- bootstrap.cf(cf$c2)
        c3 <- bootstrap.cf(cf$c3)
    }

    pdf(sprintf("output/NJN_%s_%s_%s_%s_%i_%i_to_%i.pdf",particle, g.i, g.f, mom.tag, bool, par.t1, par.t2))
    #options(warn = -1)
    print(sprintf("NJN_%s_%s_%s_%s_%i.pdf", particle, g.i, g.f, mom.tag, bool ))
    plot.cf(c2, main=paste("2pt Cor. on cA211a/b.30.32, src=", g.i, ", snk=", g.f, ", mom=", mom.tag, ",\nparticletype=", particle), log="y", ylab="C(t)", xlab="t/a", sep="" )
    plot.cf(c3, main=paste("3pt Cor. on cA211a/b.30.32, src=", g.i, ", snk=", g.f, ", mom=", mom.tag, ",\nparticletype=", particle), log="y", ylab="C(t)", xlab="t/a", sep="") #legend_title="N-J-N")
    meff <- bootstrap.effectivemass(c2, type="log")
    plotInfo$massfit <- list()
    plotInfo$massfit[[paste0("t1_",par.t1)]] <- list()
    plotInfo$masssum <- list()
    plotInfo$masssum[[paste0("t1_",par.t1)]] <- list()
    plotInfo$massfit[[paste0("t1_",par.t1)]][[paste0("t2_",par.t2)]] <- fit.effectivemass(meff, t1=par.t1, t2=par.t2) 
    plot.effectivemass(plotInfo$massfit[[paste0("t1_",par.t1)]][[paste0("t2_",par.t2)]], main=paste("Effective mass of 2pt Cor. on cA211a/b.30.32, src=", g.i, ", snk=", g.f, "\n, mom=", mom.tag, ",particletype=", particle), ylab="aM_eff", xlab="t/a", ylim=c(0, 1))
    print(plotInfo$masssum[[paste0("t1_",par.t1)]][[paste0("t2_",par.t2)]] <- summary.effectivemassfit(plotInfo$massfit[[paste0("t1_",par.t1)]][[paste0("t2_",par.t2)]]))
    #plotInfo$masssumm <- summary.effectivemassfit(meff)
    plotInfo$ratiofit <- list()
    if(is.null(plotInfo$ratiofit[[paste0("t1_",par.t1)]])){
        plotInfo$ratiofit[[paste0("t1_",par.t1)]] <- list()
    }
    if(is.null(plotInfo$ratiofit[[paste0("t1_",par.t1)]][[paste0("t2_",par.t2)]])){
        plotInfo$ratiofit[[paste0("t1_",par.t1)]][[paste0("t2_",par.t2)]] <- list()
    }
    for(t in c(1:par.tau)){
        print(paste0("Plotting comb with tau =",t))
        plotInfo$ratiofit[[paste0("t1_",par.t1)]][[paste0("t2_",par.t2)]][[paste0("tau", t)]] <- plot_ratio(c2,c3, t1=par.t1, t2=par.t2, ylim=c(-25, -5), .tau=t)
                  }

    #plot_ratiofit_vs_tau(part=particle, gi=g.i, gf=g.f, mom=num_i)
    #if(calcAll){

    #    ratioData[[particle]][[p_tot_tag]][[g.i]][[g.f]] <<- plot_ratio(c2,c3, t1=par.t1, t2=par.t2)
    #} else {
    #    ratioData <<- plot_ratio(c2,c3, t1=par.t1, t2=par.t2)
    #}
    #plot_tau(plotData, n_tau=par.tau)
    dev.off()
    return(plotInfo)
}

#' plot Function to compare 2 ratios
#'
#' This function allows you to compare two specific correlator ratios.
#' @param 
#' @keywords
#' @export
#' @examples

plot_2combs <- function(gi1 = "Gi_Cg5", gf1 = "Gf_Cg5", gi2 = "Gi_Cg5gt", gf2 = "Gf_Cg5gt", mom.tag = "px00py00pz00", particle="P", col2 = "blue"){

        stopifnot(calcAll)

        num_i <- str_count(mom.tag, "1")
        p_tot_tag <- sprintf("p_tot%i", num_i)
        c2 <- bootstrap.cf(cf_2pt[[particle]][[p_tot_tag]][[gi1]][[gf1]])
        c3 <- bootstrap.cf(cf_3pt[[particle]][[p_tot_tag]][[gi1]][[gf1]])
        
        cat("Ratio: ", gi1, " " , gf1)
        ratio1 <- plot_ratio(c2,c3, t1=par.t1, t2=par.t2)

        c2 <- bootstrap.cf(cf_2pt[[particle]][[p_tot_tag]][[gi2]][[gf2]])
        c3 <- bootstrap.cf(cf_3pt[[particle]][[p_tot_tag]][[gi2]][[gf2]])

        cat("Ratio: ", gi2, " " , gf2)
        ratio2 <- plot_ratio(c2,c3, t1=par.t1, t2=par.t2, add=TRUE, .col=col2)


        return(invisible(list("ratio1"=ratio1, "ratio2"=ratio2)))
}

#ratios <- plot_comb("Gi_Cg5","Gf_Cg5","px00py00pz00", test)
