#! /bin/bash

#######################################################################
#
#
######################################################################

###########_variables_####################
pathbase=$(realpath ./)
pathdest=$(realpath ./Simulation)
BASE_NAME=Simul
TYPE=Domestication
DATE=$(date +%F)
Compt=1
NEWSIM=$pathdest/$BASE_NAME-$TYPE-${DATE}-$Compt

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
CRS=FALSE
TIMEBOT=2800
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
	echo " -C 	 enable the switch in selection regime"
	echo " -E    allow to modify the others parameters if needed"
	echo " -m <value>    change the mutation rate Default : $MUT_RATE "
	echo " -c <value>    change the initial connection probability of the network Default : $INIT_CONNECT"
	echo " -B <value>    change the number of generations before the bottleneck Default : 12000"
	echo " -A <value>	 change the number of generations after the bottleneck  Default : 6200"
	echo " -T    allow to modify the scenario of domestication "
	echo " -D 	 register the launchfile in bigmem in order to allow the parrallel launch of all your simulation"
	echo " -h    show help"
}

	
#######_options_##########################
while getopts "h?:N:o:r:m:c:A:B:k:S::TGEDCLIZ" opt; do 
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
	C)		
	CRS=TRUE 
	OPTION=$OPTION" -crs">&2
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
	h|\?)
        show_help
        exit 0
    ;;
	\?)
	echo "option invalide, aide disposible via argument -h " >&2
	;;
esac
done
OPTION=$OPTION" -mu $MUT_RATE"
###########_actions_#######################

#creation du fichier contenant le pipeline + dossier d'entrée et sortie + sortie pour R



if [ ! -d "$pathdest" ] ; then
mkdir $pathdest
fi


while [ -d "$NEWSIM" ] ; do
Compt=$(($Compt +1 ))
NEWSIM=$pathdest/$BASE_NAME-$TYPE-${DATE}-$Compt
done

mkdir $NEWSIM
cp $pathbase/param.txt $NEWSIM
touch $NEWSIM/launchfile.sh 
chmod a+x $NEWSIM/launchfile.sh 
mkdir $NEWSIM/out-in_Fold

#modification du fichier param.txt
echo "SIMUL_OUTPUT  $OUTPUT_GEN" >> $NEWSIM/param.txt
echo "INIT_PSIZE  $SIZE_POP" >> $NEWSIM/param.txt
echo "INIT_CONNECT $INIT_CONNECT" >> $NEWSIM/param.txt
echo "GENET_MUTRATES  $MUT_RATE" >> $NEWSIM/param.txt

#creation et remplissage du fichier extparam.txt
echo "Nombre_replicat   $NB_REP" > $NEWSIM/extparam.txt
echo "Règles de codage :" >> $NEWSIM/extparam.txt
echo "1 = gène stable avec même objectif avant et après domestication" >> $NEWSIM/extparam.txt
echo "2 = gène plastique avec pour valeur x = valeur de l'environnement " >> $NEWSIM/extparam.txt
echo "3 = gène nouvellement stable (changement d'objectif quand on passe de 1 à 3 )" >> $NEWSIM/extparam.txt
echo "4 = gène plastique avec pour valeur x = 1-valeur de l'environnement" >> $NEWSIM/extparam.txt
echo "0 = gène non selectionné = gène neutre " >> $NEWSIM/extparam.txt
echo "NB : il faut autant de chiffre qu'il y a de gènes " >> $NEWSIM/extparam.txt
echo "Scenario_part1 2 1 1 1 1 1 1 2 2 2 4 4 4 0 0 0 0 0 0 0 0 0 0 0 0" >> $NEWSIM/extparam.txt
echo "$CRS"

if [ "$CRS" = "TRUE" ]; then	# if CRS = T alors le scenario 1 et le scenario 2 sont différent, sinon le scenario1=scenario2
echo "Scenario_part2 2 1 1 3 3 0 0 3 0 2 4 3 0 3 3 3 3 0 0 0 0 0 0 0 0" >> $NEWSIM/extparam.txt ; 
else  
echo "Scenario_part2 2 1 1 1 1 1 1 2 2 2 4 4 4 0 0 0 0 0 0 0 0 0 0 0 0" >> $NEWSIM/extparam.txt ; 
fi	
echo "nb generations bottleneck = $TIMEBOT" >> $NEWSIM/extparam.txt
echo "OPTIONS = $OPTION" >> $NEWSIM/extparam.txt
###
OPTION=$OPTION" -dir $NEWSIM"


# affichage et edition des fichier de paramétrage si option -T ou -E activé
cd $NEWSIM
if [ "$PARAM_MODIF" = TRUE ]; then
nano param.txt
fi
if [ "$SCENAR_MODIF" = TRUE ]; then
nano extparam.txt 
fi

#lancement de la géneration des fichier d'entrée + création du launchfile
echo "création de l'environnement"
cd $pathbase/Script
echo "$OPTION"
./SimFiles_generator.R $OPTION
#register the lauchfile.sh of newsim in bigmem.sh
if [ "$BIG_LAUNCH" = TRUE ]; then
cat "$NEWSIM/launchfile.sh" >> $pathbase/bigmem.sh
echo " le fichier de lancement de simulation à été transféré dans bigmem"
fi

echo "fin de tache"
