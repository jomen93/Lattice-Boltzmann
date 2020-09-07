
clear



# define directories
Animation_folder=Animations
Data_folder=DATA
Plots_folder=plots

# verification of parents folders 
if [ ! -e "${Animation_folder}" ]; then
	mkdir "${Animation_folder}"
	echo "The folder ${Animation_folder} does not exist, it's created" 
else
	rm "${Animation_folder}"/*.gif
	echo "Old elements in ${Animation_folder} are deleted"
fi

if [ ! -e "${Data_folder}" ]; then
	mkdir "${Data_folder}"
	echo "The folder ${Data_folder} does not exist, it's created" 
else
	rm "${Data_folder}"/*.dat
	echo "Old elements in ${Data_folder} are deleted"
fi

if [ ! -e "${Plots_folder}" ]; then
	mkdir "${Plots_folder}"
	echo "The folder ${Plots_folder} does not exist, it's created" 
else
	rm "${Plots_folder}"/*.png
	echo "Old elements in ${Plots_folder} are deleted"
fi


#compilation of principal program
python src/main.py

# organization of utput files
mv *.gif "${Animation_folder}"
mv *.dat "${Data_folder}"
mv *.png "${Plots_folder}"







