library(nGn.Analysis)


gi= "Gi_Cg5"
gf = "Gf_Cg5"
#gi= "Gi_Cg5gt"
#gf = "Gf_Cg5gt"
mom_tag <- "px00py00pz00"
max_conf <- 1756   
min_conf <- 4
stepwidth <- 4
test <- FALSE
calcAll <- TRUE
doGEVP <- FALSE
doTau <- FALSE
par.tau <- 8
time <- 64
par.t1 <- 4
par.t2 <- 12

if(test){
    max_conf <-100
    min_conf <- 4
    stepwidth <- 24
}

num_conf <- (max_conf-min_conf)/stepwidth+1

if(calcAll){

    mom_tag <- list(   # "px01py01pz-01",
                       # "px01py01pz01",
                       # "px01py-01pz-01",
                       # "px01py-01pz01",
                       # "px-01py01pz-01",
                       # "px-01py01pz01",
                       # "px-01py-01pz-01",
                       # "px-01py-01pz01",
                       # "px01py01pz00",
                       # "px01py-01pz00",
                       # "px-01py01pz00",
                       # "px-01py-01pz00",
                       # "px01py00pz-01",
                       # "px01py00pz01",
                       # "px-01py00pz-01",
                       # "px-01py00pz01",
                       # "px00py01pz-01",
                       # "px00py01pz01",
                       # "px00py-01pz-01",
                       # "px00py-01pz01",
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
                    #c( "Gi_Cgt","Gf_Cgt"),
                    #c( "Gi_C","Gf_C"),
                    c( "Gi_Cg5","Gf_Cg5gt"),
                    c( "Gi_Cg5gt","Gf_Cg5")
                    
                                            )
        cf_2pt <- list()
        cf_3pt <- list()
        
        calc.all(test, "ab")
        #calc.all(test, "b")

} else {
    #d[[sprintf("%s_%s_%s", gi, gf, mom_tag)]] <- calc.cf(gi, gf, mom_tag, test)
            cf <- calc.cf(gi, gf, mom_tag, test)
            
}

ratios <- plot.comb("Gi_Cg5","Gf_Cg5","px00py00pz00", test)
