all :
	gcc -ggdb -Wall -o fd-solver domain1d.c main.c node1d.c update-functions.c material-data/can/can.c material-data/can/datafile.c -Imaterial-data/can -lm

