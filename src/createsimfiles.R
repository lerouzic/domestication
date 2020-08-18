#!/usr/bin/env Rscript

baseloc<-getwd() 


###_defualt_values_###
G1<-12000
G2<-6200
k<-2.45
bottleneck_time<-2800
digits <- 4
bottleswitch<-"OFF"
restrict_env<-"OFF"
selectionswitch<-"OFF"
standar_dev<-"OFF"
strengthswitch<-"OFF"
delta_env<-0.2
zero_selec<-"OFF"
###_options_####
arg_randopt<-commandArgs(trailingOnly = TRUE)
arg_count<-0
for (arg in arg_randopt){
arg_count<-arg_count+1
if(arg=="-dir"){simudir<-arg_randopt[arg_count+1]}
if (arg=="-h") {print(domest_help) ; suppressMessages(stop())}
if (arg=="-bot"){bottleswitch<-"ON"}
if (arg=="-ren"){restrict_env<-"ON"}
if (arg=="-crs"){selectionswitch<-"ON"}
if (arg=="-G1"){G1<-as.numeric(arg_randopt[arg_count+1])}
if (arg=="-G2"){G2<-as.numeric(arg_randopt[arg_count+1])}
if (arg=="-k"){k<-arg_randopt[arg_count+1]}
if (arg=="-env"){restrict_env<-"ON"}
if (arg=="-mu"){init_mutrate<-as.numeric(arg_randopt[arg_count+1])}
if (arg=="-std"){standar_dev<-"ON"}
if (arg=="-strengthswitch"){strengthswitch<-"ON";newstrength<-as.numeric(arg_randopt[arg_count+1])}
if (arg=="-zero_selec"){zero_selec<-"ON"}
}
if (dir.exists(normalizePath(simudir))==T){simudir<-normalizePath(simudir)} else {stop("The file doesn't exist")}

####_sources_###
source ("sampleur.R")
source ("readfiles.R")
source("Function_repository.R")

###_destinations_###
setwd(simudir)
z<-simudir
param <- "param.txt"
launch <- "./launchfile.sh"
simulpath <- "~/SimuDom/Simul_Prog"  

##########_help_text_#############
domest_help<-c(" -h show the help",
             "-bo change the population size over time during domestication",
             "-ren retrain the environement fluctuation during domestication")
####################
load("varsave.RData")
bottleneck_time<-genbot

print(paste(bottleswitch,selectionswitch,G1,G2,k))
bottleneck_popsize<-bottleneck_time*as.numeric(k)/2
pop_vector<-rep(popsizemax,3)
if (bottleswitch=="ON") {
pop_vector<-c(popsizemax,bottleneck_popsize,popsizemax)
}

time_vector<-c(G1,bottleneck_time,G2)

simul_maxgen<-sum(time_vector)
cat(paste0("SIMUL_MAXGEN\t",simul_maxgen),file="param.txt",append=T,sep = "\n")


#####_parameter_creating_process_start_########

for ( a in 1:rep) {
dir_rep_file<-paste0(simudir,"/out-in_Fold","/","rep_file",a)
dir.create(dir_rep_file)


ss<-dir_rep_file
sss <- paste0(z,"/out-in_Fold")



#######_variables_########
optimumtext<-vector()
popsizetext<-vector()
strenghttext<-vector()
scenariotext<-vector()
mutratetext<-vector()

#######

for (case in 1:length(time_vector)) {
  popsizetext<-c(popsizetext,rep(paste(pop_vector[case],collapse = "\t"),time_vector[case]))
}



#######
optmem_part1<-sapply(domesticatio_scenario_part1 , opt_creator,type=1 )
optmem_part2<-optmem_part1*sapply(domesticatio_scenario_part2,switchoff_opt)+sapply( domesticatio_scenario_part2, opt_creator,type=3)

vec_scenario<-domesticatio_scenario_part1
stable_opt<-optmem_part1
part_count<-0
if (zero_selec=="OFF"){strength_of_selection<-10} else {strength_of_selection<-0}


for (part in 1:simul_maxgen) {
  
  part_count<-part_count +1
  
  if (part_count==(G1+1)){
    vec_scenario<-domesticatio_scenario_part2
    stable_opt<-optmem_part2
    if (strengthswitch=="ON"){strength_of_selection<-newstrength}
  }
  if (standar_dev =="ON" & part_count > G1){mutrate<-0} else {mutrate<-init_mutrate}
  env<-round(runif(1),digits = 4)
  if (restrict_env=="ON" & part_count > G1){env<-round(runif(1,0.5-delta_env,0.5+delta_env),digits = 4)}
  mutratetext[length(mutratetext)+1]<-mutrate
  strenghttext[length(strenghttext)+1]<-paste(sapply(vec_scenario, strength_switch,strength_of_selection=strength_of_selection) , collapse = "\t")
  optimumtext[length(optimumtext)+1]<-paste(stable_opt+sapply( vec_scenario, plastic_optcreator,env=env) ,collapse = "\t")
  scenariotext[length(scenariotext)+1]<-paste(vec_scenario,collapse = "\t")
}


########_creating_the_files_##########
writeLines(optimumtext,paste0(ss,"/optimum.txt"))
writeLines(popsizetext,paste0(ss,"/pop_size.txt"))
writeLines(strenghttext,paste0(ss,"/strenght_select.txt"))
writeLines(scenariotext,paste0(ss,"/scenario_per_rep.txt"))
writeLines(as.character(mutratetext),paste0(ss,"/mutrate.txt"))
optgenerator(opt=paste0(ss, "/optimum.txt"),popsize=paste0(ss,"/pop_size.txt"),select_strenght=paste0(ss,"/strenght_select.txt"),mutrate=paste0(ss,"/mutrate.txt"), param=param, i=a,sss)


#######_creating_the_launcher_########

cat(simulpath, " -p ", sss, "/par_rep", a , "_gen1 -o ", sss, "/rep", a , ".txt\n", sep="", file=launch, append=T)



}

