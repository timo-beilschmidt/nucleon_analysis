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
recalc <- TRUE
test <- FALSE
calcAll <- TRUE
doGEVP <- FALSE
doTau <- FALSE
doPlot <- FALSE
par.tau <- 8
time <- 64
par.t1 <- 6
par.t2 <- 12
ratios <- list()

if(test){
    max_conf <-100
    min_conf <- 4
    stepwidth <- 24
}

num_conf <- (max_conf-min_conf)/stepwidth+1
    if(calcAll){

        mom_tag <- list(    "px01py01pz-01",
                            "px01py01pz01",
                            "px01py-01pz-01",
                            "px01py-01pz01",
                            "px-01py01pz-01",
                            "px-01py01pz01",
                            "px-01py-01pz-01",
                            "px-01py-01pz01",
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
                        c( "Gi_C","Gf_C"),
                        c( "Gi_Cg5","Gf_Cg5gt"),
                        c( "Gi_Cg5gt","Gf_Cg5")
                        
                                                )
        if(recalc){

            cf_2pt <- list()
            cf_3pt <- list()
            ratioData <- list() 
            calc_all(test, "ab")
            #calc_all(test, "b")
            sink("./log.txt")
            for(p_tag in c(0:1)){
                p_tot_tag <- sprintf("p_tot%i", p_tag)
                if(p_tag>0){
                    num <- strrep("1", p_tag)
                } else {
                    num <- "0"
                }
                if(is.null(ratios[[p_tot_tag]])){
                    ratios[[p_tot_tag]] <- list()

                }
                for( g in gi){
                    if(is.null(ratios[[p_tot_tag]][[g[1]]])){
                        ratios[[p_tot_tag]][[g[1]]] <- list()
                    }
                    ratios[[p_tot_tag]][[g[1]]][[g[2]]] <- plot_comb(g[1], g[2], num, test ) 

                }
            }
            sink()

        } 
    
    } else {
        #d[[sprintf("%s_%s_%s", gi, gf, mom_tag)]] <- calc_cf(gi, gf, mom_tag, test)
                cf <- calc_cf(gi, gf, mom_tag, test)
                
    }

if(doPlot){
    
    #ratios$a <- plot_comb("Gi_Cg5","Gf_Cg5","px00py00pz00", test)
    #ratios$b <- plot_comb("Gi_Cg5gt","Gf_Cg5gt","px00py00pz00", test)
    #ratios$c <- plot_comb("Gi_Cg5","Gf_Cg5gt","px00py00pz00", test)
    #ratios$d <- plot_comb("Gi_Cg5gt","Gf_Cg5","px00py00pz00", test)
    #ratios$e <- plot_comb("Gi_Cg5","Gf_Cg5","px01py00pz00", test)
    #ratios$f <- plot_comb("Gi_Cg5gt","Gf_Cg5gt","px01py00pz00", test)
    #ratios$g <- plot_comb("Gi_Cg5","Gf_Cg5gt","px01py00pz00", test)
    #ratios$h <- plot_comb("Gi_Cg5gt","Gf_Cg5","px01py00pz00", test)
    #ratios$aa <- plot_comb("Gi_C","Gf_C","px00py00pz00", test)
    #ratios$bb <- plot_comb("Gi_Cgt","Gf_Cgt","px00py00pz00", test)

    plot_vs <- function(gi1 = "Gi_Cg5", gf1 = "Gf_Cg5gt", gi2 = "Gi_Cg5gt", gf2 = "Gf_Cg5", mom.tag = "p0", col2 = "blue"){

        pdf(sprintf("output/NJN_%s_%s_vs_%s_%s_%s_%i.pdf", gi1, gf1, gi2, gf2, mom.tag, test))
        r <- plot_2combs(gi1, gf1, gi2, gf2, mom.tag, col2)
        dev.off()
        return(r)
    }


    #ratios$i <- plot_vs(gi1 = "Gi_Cg5", gf1 = "Gf_Cg5", gi2 = "Gi_Cg5gt", gf2 = "Gf_Cg5gt", mom.tag = "p0", col2 = "blue")

    #ratios$j <- plot_vs(gi1 = "Gi_Cg5", gf1 = "Gf_Cg5", gi2 = "Gi_Cg5gt", gf2 = "Gf_Cg5gt", mom.tag = "p1", col2 = "blue")

    #ratios$k <- plot_vs(gi1 = "Gi_Cg5", gf1 = "Gf_Cg5gt", gi2 = "Gi_Cg5gt", gf2 = "Gf_Cg5", mom.tag = "p0", col2 = "blue")

    #ratios$l <- plot_vs(gi1 = "Gi_Cg5", gf1 = "Gf_Cg5gt", gi2 = "Gi_Cg5gt", gf2 = "Gf_Cg5", mom.tag = "p1", col2 = "blue")

}
