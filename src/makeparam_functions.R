Rutils <- require(R.utils, quietly=TRUE) # Path manipulation 

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

