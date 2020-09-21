# To be "Sourced" from a R script file

Rutils <- suppressPackageStartupMessages(require(R.utils)) # Path manipulation 

create.launchfile <- function(prog.path, param.files, output.files, launch.file="./launchfile.sh", relative.paths=TRUE) {
	# relative.paths sets all paths relative to the launch.file directory
	stopifnot(
		length(param.files) > 0, 
		length(param.files) == length(output.files),
		file.exists(prog.path), 
		all(file.exists(param.files)))
		
	if (Rutils && relative.paths) {
		launch.dir <- dirname(launch.file)
		prog.path   <- getRelativePath(prog.path, launch.dir)
		param.files <- getRelativePath(param.files, launch.dir)
		output.files<- getRelativePath(output.files, launch.dir) 
	}
	
	commands <- paste0(prog.path, " -p ", param.files, " -o ", output.files)
	writeLines(commands, launch.file)
}

read.param <- function(param.file) {
	# The structure of the parameter file is assumed to be, for each line:
	# PARNAME x y
	# The function returns a list named by the parameter names
	
	lines <- readLines(param.file)
	lines <- lines[!grepl(lines, pattern="^#")]
	splines <- strsplit(lines, split="\\s+")
	nn <- sapply(splines, "[", 1)
	ans <- lapply(splines, function(ll) {
		if (length(ll) <= 1) return(NA)
		ll <- ll[-1]
		numll <- suppressWarnings(as.numeric(ll))
		if (any(is.na(numll))) 
			return(ll)
		else 
			return(numll)
		})
	names(ans) <- nn
	ans
}

write.param <- function(param.file, param) {
	content <- paste(names(param), sapply(param, paste, collapse="\t"), sep="\t")
	writeLines(content, param.file)
}

make.randopt <- function(oldopt, pattern, begin.bottleneck=FALSE) {
	# Returns a vector of random optima according to the provided pattern
	# (constant vs fluctuating optima)
	stopifnot(length(oldopt) == length(pattern),
			all(pattern %in% 0:5))
	newopt <- oldopt
	env <- runif(1)
	newopt[pattern == 0] <- 0
	newopt[pattern == 1] <- oldopt[pattern == 1] # does nothing, but clearer this way
	newopt[pattern == 2] <- env
	if (begin.bottleneck) newopt[pattern == 3] <- runif(sum(pattern == 3))
	newopt[pattern == 4] <- 1-env
	newopt[1] <- env
	
	newopt
}

make.selstr <- function(pattern) {
	# Mirrors the make.random function, but for the selection strength. 
	# The algorithm is trivial, just sets zero for non-selected genes
	# (including the first one = signal). 
	
	ans <- pattern != 0
	ans[1] <- 0
	ans
}

create.paramseries <- function(param.template.file, extparam.file, simul.dir, overwrite=FALSE, verbose=FALSE, allow.extrapar=c("GENET_MUTRATES")) 
	# This is the main algorithm that create the simulation structure. 
	# The function retruns the necessary information to make a launchfile
	# (it does not write the launchfile, because all information is not available here)
	# Two parameter files are necessary : 
	#  - the parameter template, which is a regular parameter file that will be used for the default settings
	#  - the "extended" parameter file, which provides "meta" information to run the simulation program
	#    this is specific for the domestication project, and involves three simulation stages:
	#     . before the bottleneck (burn-in)
	#     . during the bottleneck
	#     . after the bottleneck
	# It is expected that the information from the two parameter files match (e.g. number of genes), 
	# consequences of a non-matching pattern are undefined. 
	# The algorithm reads the two parameter files and create the necessary files for simulations
	# For each simulation replicate (REPLICATES variable in the ext.par file), a simulation directory
	# is created, and as many parameter files as generations are created within each directory
	#  Note: this represents a very large number of (small) files, which can represent a burder for 
	#  the system, especially when creating the files or compressing the directory
	# The first parameter file of a simulation is complete, the next ones only change what is necessary:
	#   . every generation: the environmental variable and the fluctuating optima
	#   . before and after the bottleneck: the population size
	#   . at the beginning of the bottleneck: the "domestication" selection pattern
	# The directory and parameter naming pattern is specified in a few internal functions. 

	 {
		
	.repID <- function(rep, ndigits=4)
		formatC(rep, width = ndigits, format = "d", flag = "0")
		
	.genID <- function(gen, ndigits=5)
		formatC(gen, width = ndigits, format = "d", flag = "0")
		
	.repDir  <- function(rep, rep.template = "rep") 
		paste0(rep.template, "-", .repID(rep))
		
	.repFile <- function(rep, gen, rep.template = "rep", gen.template="gen")
		paste0("param-", rep.template, .repID(rep), "-", gen.template, .genID(gen), ".txt")
		
	.outFile <- function(rep, rep.template = "rep")
		paste0("out-", rep.template, .repID(rep), ".txt")
	
	stopifnot(file.exists(param.template.file), file.exists(extparam.file))

	param.template <- read.param(param.template.file)
	extparam       <- read.param(extparam.file)
	
	myparam <- param.template
	
	extrapar <- extparam[allowed.extrapar]
	
	# A few shortcuts to make the code more readable
	bo.b <- extparam$BOTTLENECK_BEGIN 
	bo.d <- extparam$BOTTLENECK_DURATION
	bo.a <- extparam$BOTTLENECK_AFTER
	totgen <- bo.b + bo.d + bo.a
	sel.strength <- max(myparam$FITNESS_STRENGTH)
	sel.strength2 <- if ("FITNESS_STRENGTH" %in% names(extparam)) max(extparam$FITNESS_STRENGTH) else sel.strength
	
	# Checks for the consistency of both parameter files
	if (!"SIMUL_MAXGEN" %in% names(myparam)) 
		myparam$SIMUL_MAXGEN <- totgen
		
	if (myparam$SIMUL_MAXGEN < totgen) 
		warning("Inconsistency between the MAXGEN parameters in ", param.template.file, 
		        " and bottleneck duration in ", extparam.file, ".")
	
	
	stopifnot(
		length(extparam$SCENARIO_PART1) == myparam$GENET_NBLOC , 
		length(extparam$SCENARIO_PART1) == length(extparam$SCENARIO_PART2))
	
	stopifnot(extparam$REPLICATES > 0)
	
	# Now the parameter files:
	dir.create(file.path(simul.dir), showWarnings = FALSE)
	
	# It looks difficult to make it simpler: parameter files need to be updated
	# on a case-by-case basis, as the parameters to change depend on the 
	# simulation stage (before, during, of after the bottleneck)
	
	if (verbose)
		pb <- txtProgressBar(min = 1, max = extparam$REPLICATES, style = 3)
	for (rep in 1:extparam$REPLICATES) {
		if (verbose) setTxtProgressBar(pb, rep)
		
		dir.create(file.path(simul.dir, .repDir(rep)), showWarnings = FALSE)
		if (overwrite) { # This is quite powerful, use with caution (simulation results are deleted prior to launching a new sim)
			unlink(list.files(path=file.path(simul.dir, .repDir(rep)), full.names=TRUE))
		}
		
		# First generation: use the full template file
		optim <- make.randopt(runif(myparam$GENET_NBLOC), extparam$SCENARIO_PART1)
		myparam$FITNESS_OPTIMUM <- optim
		myparam$FITNESS_STRENGTH <- sel.strength*make.selstr(extparam$SCENARIO_PART1)
		myparam$FILE_NEXTPAR <- suppressWarnings(normalizePath(file.path(simul.dir, .repDir(rep), .repFile(rep, 1))))
		par.file.name <- file.path(simul.dir, .repDir(rep), .repFile(rep, 0))
		if (!file.exists(par.file.name) || overwrite)
			write.param(par.file.name, myparam)
		
		# From here, just update what is necessary
		if (bo.b > 2)
		for (gen in 1:(bo.b-1)) {
			optim <- make.randopt(optim, extparam$SCENARIO_PART1)
			par.file.name <- file.path(simul.dir, .repDir(rep), .repFile(rep, gen))
			if (!file.exists(par.file.name) || overwrite)
				write.param(par.file.name,
					list(FITNESS_OPTIMUM = optim, 
						 FILE_NEXTPAR    = suppressWarnings(normalizePath(file.path(simul.dir, .repDir(rep), .repFile(rep, gen+1))))))
		}

		# first bottleneck generation
		optim <- make.randopt(optim, extparam$SCENARIO_PART2, begin.bottleneck=TRUE)
		par.file.name <- file.path(simul.dir, .repDir(rep), .repFile(rep, bo.b))
		if (!file.exists(par.file.name) || overwrite)
			write.param(par.file.name,
					c(list(INIT_PSIZE = round(bo.d*extparam$BOTTLENECK_STRENGTH/2), # Maize is diploid, so the strenght k = 2N/t
						FITNESS_OPTIMUM = optim, 
						FITNESS_STRENGTH     = sel.strength2*make.selstr(extparam$SCENARIO_PART2),
						FILE_NEXTPAR    = suppressWarnings(normalizePath(file.path(simul.dir, .repDir(rep), .repFile(rep, bo.b+1))))),
						extrapar))
					 
		# The rest of the bottleneck
		if (bo.d > 1)
		for (gen in ((bo.b+1):(bo.b+bo.d))) {
			optim <- make.randopt(optim, extparam$SCENARIO_PART2)
			par.file.name <- file.path(simul.dir, .repDir(rep), .repFile(rep, gen))
			if (!file.exists(par.file.name) || overwrite)
				write.param(par.file.name,
					list(FITNESS_OPTIMUM = optim, 
						FILE_NEXTPAR    = suppressWarnings(normalizePath(file.path(simul.dir, .repDir(rep), .repFile(rep, gen+1))))))
		}
		
		# First generation after the bottleneck
		if (bo.a > 1){
		optim <- make.randopt(optim, extparam$SCENARIO_PART2)
		par.file.name <- file.path(simul.dir, .repDir(rep), .repFile(rep, bo.b + bo.d))
		if (!file.exists(par.file.name) || overwrite)
			write.param(par.file.name,
				list(FITNESS_OPTIMUM = optim, 
					INIT_PSIZE      = myparam$INIT_PSIZE,
					FILE_NEXTPAR    = suppressWarnings(normalizePath(file.path(simul.dir, .repDir(rep), .repFile(rep, bo.b + bo.d + 1))))))
		}
		
		# Domestication after bottleneck
		if (bo.a > 2)
		for (gen in ((bo.b+bo.d+1):(bo.b+bo.d+bo.a-1))) {
			optim <- make.randopt(optim, extparam$SCENARIO_PART2)
			par.file.name <- file.path(simul.dir, .repDir(rep), .repFile(rep, gen))
			if (!file.exists(par.file.name) || overwrite)
				write.param(file.path(simul.dir, .repDir(rep), .repFile(rep, gen)),
					list(FITNESS_OPTIMUM = optim, 
						FILE_NEXTPAR    = suppressWarnings(normalizePath(file.path(simul.dir, .repDir(rep), .repFile(rep, gen+1))))))
		}
		
		# Very last generation (no NEXTPAR)
		optim <- make.randopt(optim, extparam$SCENARIO_PART2)
		
		write.param(file.path(simul.dir, .repDir(rep), .repFile(rep, bo.b+bo.d+bo.a)),
				list(FITNESS_OPTIMUM = optim))
	}
	
	# Returns the parameter (and output) file names that will be necessary to make the launchfile
	return(list(param=file.path(simul.dir, .repDir(1:extparam$REPLICATES), .repFile(1:extparam$REPLICATES, 0)), 
	            out =file.path(simul.dir, .repDir(1:extparam$REPLICATES), .outFile(1:extparam$REPLICATES))))
}
