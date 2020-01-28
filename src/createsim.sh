#! /bin/bash

#######################################################################
#
#
######################################################################

###########_PATH_#########################
DIR_ROOT=../
DIR_SRC=$DIR_ROOT/src/
DIR_SIM=$DIR_ROOT/cache/
DIR_PAR=$DIR_ROOT/param/
PATH_GENERATOR=$DIR_SRC/SimFiles_generator.R

###########_variables_####################
BASE_NAME=Simul
TYPE=Domestication
DATE=$(date +%F)
Compt=1
NEWSIM=$DIR_SIM/$BASE_NAME-$TYPE-${DATE}-$Compt

#default values
SIZE_POP=1000
OUTPUT_GEN=100
NB_REP=8
NB_PLAST_GEN=6
NB_SELECT_GEN=12
MUT_RATE=0.001
INIT_CONNECT=1
BIG_LAUNCH=FALSE 
DOMEST=FALSE
SCENAR_MODIF=FALSE
PARAM_MODIF=FALSE
OPTION=I
TIMEBOT=2800
TEMPLATE_PARAM=$DIR_PAR/param0.txt
TEMPLATE_EXTPARAM=$DIR_PAR/extparam0.txt

#########_thelp_function_######################
show_help() {
	echo " -Z    switch the strength of selection to 0 during all the time of the simulation"
	echo " -I    switch the mutation rate to 0 at the begining of the domestication"
	echo " -S <value>	 switch the strength of the selection to <value> at the begining of the domestication part "
	echo " -L 	 reduce the environnement varaition at the begining of the domestication part"
	echo " -N <value>    allow to modify the population size : $SIZE_POP "
	echo " -G    enable the genetic bottleneck "
	echo " -o <value>    modify the output rate  Default : $OUTPUT_GEN "
	echo " -r <value>    modify the number of replicates  Default : $NB_REP "
	echo " -k <value>	 modify the stregth of the bottleneck "
	echo " -m <value>    change the mutation rate Default : $MUT_RATE "
	echo " -c <value>    change the initial connection probability of the network Default : $INIT_CONNECT"
	echo " -B <value>    change the number of generations before the bottleneck Default : 12000"
	echo " -A <value>	 change the number of generations after the bottleneck  Default : 6200"
	echo " -D 	 register the launchfile in bigmem in order to allow the parrallel launch of all your simulation"
	echo " -P <value>    template parameter file"
	echo " -p <value>    template extended parameter file"
	echo " -h    show help"
}

	
#######_options_##########################
while getopts "h?:N:o:r:m:c:A:B:k:S::TGEDCLIZP:p:" opt; do 
case $opt in
	N)		
	SIZE_POP=$OPTARG >&2
	;;
	A)		
	OPTION=$OPTION" -G2 $OPTARG" >&2
	;;
	B)		
	OPTION=$OPTION" -G1 $OPTARG" >&2
	;;
	k)		
	OPTION=$OPTION" -k $OPTARG" >&2
	;;
	G)		
	OPTION=$OPTION" -bot" >&2
	;;
	o)		
	OUTPUT_GEN=$OPTARG >&2
	;;
	r) 		
	NB_REP=$OPTARG >&2
	;;
	m) 		
	MUT_RATE=$OPTARG >&2
	;;
	c) 		
	INIT_CONNECT=$OPTARG >&2
	;;
	T) 		
	SCENAR_MODIF=TRUE >&2
	;;
	E)		
	PARAM_MODIF=TRUE >&2
	;;
	D)		
	BIG_LAUNCH=TRUE >&2
	;;
	L)		
	OPTION=$OPTION" -env" >&2
	;;
	S)		
	OPTION=$OPTION" -strengthswitch $OPTARG" >&2
	;;
	I)		
	OPTION=$OPTION" -std" >&2
	;;
	Z)      
	OPTION=$OPTION" -zero_selec" >&2
	;;
	P)
	TEMPLATE_PARAM=$OPTARG >&2
	;;
	p)
	TEMPLATE_EXTPARAM=$OPTARG >&2
	;;
	h|\?)
        show_help
        exit 0
    ;;
	\?)
	echo "Invalid option, help avilable via the argument -h" >&2
	;;
esac
done
OPTION=$OPTION" -mu $MUT_RATE"


###########_actions_#######################

# Creating the directory including the simulation pipeline + input/output files


if [ ! -d "$DIRSIM" ] ; then
mkdir $DIRSIM
fi


while [ -d "$NEWSIM" ] ; do
Compt=$(($Compt +1 ))
NEWSIM=$DIRSIM/$BASE_NAME-$TYPE-${DATE}-$Compt
done

SIMPARAM=$NEWSIM/param.txt
SIMEXTPARAM=$NEWSIM/extparam.txt

mkdir $NEWSIM
cp $TEMPLATE_PARAM $SIMPARAM
cp $TEMPLATE_EXTPARAM $SIMEXTPARAM
touch $NEWSIM/launchfile.sh 
chmod a+x $NEWSIM/launchfile.sh 
mkdir $NEWSIM/out-in_Fold

# Changes in the parameter file
echo "SIMUL_OUTPUT  $OUTPUT_GEN"  >> $SIMPARAM
echo "INIT_PSIZE  $SIZE_POP"      >> $SIMPARAM
echo "INIT_CONNECT $INIT_CONNECT" >> $SIMPARAM
echo "GENET_MUTRATES  $MUT_RATE"  >> $SIMPARAM

# Changes in the extended parameter file

echo "Nombre_replicat   $NB_REP"            >> $SIMEXTPARAM
echo "nb generations bottleneck = $TIMEBOT" >> $SIMEXTPARAM
echo "OPTIONS = $OPTION"                    >> $SIMEXTPARAM

###
OPTION=$OPTION" -dir $NEWSIM"


# Creating the simulation environment
#~ cd $pathbase/Script
Rscript $PATH_GENERATOR $OPTION # TODO: output in the right directory
#register the lauchfile.sh of newsim in bigmem.sh
#~ if [ "$BIG_LAUNCH" = TRUE ]; then
#~ cat "$NEWSIM/launchfile.sh" >> $pathbase/bigmem.sh
#~ echo " le fichier de lancement de simulation à été transféré dans bigmem"
#~ fi
