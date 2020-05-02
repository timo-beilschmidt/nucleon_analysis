library(hadron)
library(Raff)
library(boot)
library(stringr)
library(pracma)
library(DescTools)
library(binhf)
library(nGn.Analysis)
source('/qbigwork2/beilschmidt/code/R/get_sl.R')
source('/qbigwork2/beilschmidt/code/R/R_gamma.R')

gi= "Gi_Cg5"
gf = "Gf_Cg5"
#gi= "Gi_Cg5gt"
#gf = "Gf_Cg5gt"
mom_tag <- "px00py00pz00"
max_conf <- 1756   
min_conf <- 4
stepwidth <- 4
recalc <- FALSE
test <- FALSE
calcAll <- TRUE
doGEVP <- FALSE
doTau <- FALSE
doPlot <- TRUE
par.tau <- 8
time <- 64
par.t1 <- 6
par.t2 <- 12

if(test){
    max_conf <-100
    min_conf <- 4
    stepwidth <- 24
}

#num_conf <- (max_conf-min_conf)/stepwidth+1
    if(calcAll){

        mom_tag <- list(   # "px01py01pz-01",
                           # "px01py01pz01",
                           # "px01py-01pz-01",
                           # "px01py-01pz01",
                           # "px-01py01pz-01",
                           # "px-01py01pz01",
                           # "px-01py-01pz-01",
                           # "px-01py-01pz01",
                            "px01py01pz00",
                            "px01py-01pz00",
                            "px-01py01pz00",
                            "px-01py-01pz00",
                            "px01py00pz-01",
                            "px01py00pz01",
                            "px-01py00pz-01",
                            "px-01py00pz01",
                            "px00py01pz-01",
                            "px00py01pz01",
                            "px00py-01pz-01",
                            "px00py-01pz01",
                            "px01py00pz00",
                            "px-01py00pz00",
                            "px00py01pz00",
                            "px00py-01pz00",
                            "px00py00pz-01",
                            "px00py00pz01",
                            "px00py00pz00"
                                                  )
        gi <- list(
                        c("Gi_Cg5","Gf_Cg5"),
                        c( "Gi_Cg5gt","Gf_Cg5gt"),
                        c( "Gi_Cgt","Gf_Cgt"),
                       # c( "Gi_C","Gf_C"),
                       # c( "Gi_Cg5","Gf_Cg5gt"),
                       # c( "Gi_Cg5gt","Gf_Cg5"),
                       # c( "Gi_C","Gf_Cgt"),
                       # c( "Gi_Cgt","Gf_C")

                                                )
        if(recalc){

            #ratios <- list()
            #cf_2pt <- list()
            #cf_3pt <- list()
            #ratioData <- list() 
            calc_all(test, "ab", particle="P")
            calc_all(test, "b", particle="N")
            #calc_all(test, "b")
           # sink("./log.txt")
         #   for(part in c("P", "N")){

         #       for(p_tag in c(0:1)){

         #           p_tot_tag <- sprintf("p_tot%i", p_tag)
         #           if(p_tag>0){
         #               num <- strrep("1", p_tag)
         #           } else {
         #               num <- "0"
         #           }
         #           
         #           if(is.null(ratios[[part]])){
         #               ratios[[part]] <- list()

         #           }

         #           if(is.null(ratios[[part]][[p_tot_tag]])){
         #               ratios[[part]][[p_tot_tag]] <- list()

         #           }
         #           for( g in gi){
         #               if(is.null(ratios[[part]][[p_tot_tag]][[g[1]]])){
         #                   ratios[[part]][[p_tot_tag]][[g[1]]] <- list()
         #               }
         #               ratios[[part]][[p_tot_tag]][[g[1]]][[g[2]]] <- plot_comb(g[1], g[2], num, test, particle=part ) 

         #           }
         #       }
         #   #    sink()

         #   } 
        }
    
    } else {
        #d[[sprintf("%s_%s_%s", gi, gf, mom_tag)]] <- calc_cf(gi, gf, mom_tag, test)
                cf <- calc_cf(gi, gf, mom_tag, test)
                
    }

if(doPlot){
    

    plot_vs <- function(gi1 = "Gi_Cg5", gf1 = "Gf_Cg5gt", gi2 = "Gi_Cg5gt", gf2 = "Gf_Cg5", mom.tag = "p0", col2 = "blue"){

        pdf(sprintf("output/NJN_%s_%s_vs_%s_%s_%s_%i.pdf", gi1, gf1, gi2, gf2, mom.tag, test))
        r <- plot_2combs(gi1, gf1, gi2, gf2, mom.tag, col2)
        dev.off()
        return(r)
    }
    

    for(t1 in c(4:9)){
        snd <- c((t1+3):15)
        par.t1 <- t1
        for(t2 in snd){
            par.t2 <- t2
            sink(paste0("./log",par.t1, "_to_", par.t2 ,".txt"))

    #ratios$i <- plot_vs(gi1 = "Gi_Cg5", gf1 = "Gf_Cg5", gi2 = "Gi_Cg5gt", gf2 = "Gf_Cg5gt", mom.tag = "p0", col2 = "blue")

    #ratios$j <- plot_vs(gi1 = "Gi_Cg5", gf1 = "Gf_Cg5", gi2 = "Gi_Cg5gt", gf2 = "Gf_Cg5gt", mom.tag = "p1", col2 = "blue")

    #ratios$k <- plot_vs(gi1 = "Gi_Cg5", gf1 = "Gf_Cg5gt", gi2 = "Gi_Cg5gt", gf2 = "Gf_Cg5", mom.tag = "p0", col2 = "blue")

    #ratios$l <- plot_vs(gi1 = "Gi_Cg5", gf1 = "Gf_Cg5gt", gi2 = "Gi_Cg5gt", gf2 = "Gf_Cg5", mom.tag = "p1", col2 = "blue")
            for(part in c("P", "N")){

                for(p_tag in c(0:2)){

                    p_tot_tag <- sprintf("p_tot%i", p_tag)
                    if(p_tag>0){
                        num <- strrep("1", p_tag)
                    } else {
                        num <- "0"
                    }
                    for( g in gi){
                        ratios[[part]][[p_tot_tag]][[g[1]]][[g[2]]] <- plot_comb(g[1], g[2], num, test, particle=part ) 
                        #nGn.Analysis:::plot_ratiofit_vs_tau(part=part, .gi=g[1], .gf=g[2], mom=str_count(num, "1"))
                    }
                }
            }

                sink()
        }
    }
}


