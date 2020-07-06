
data=/home/g/gabriele/pfs/targets/

while read l; do 

echo "##### $l #####"

tmpfolder=/scratch/T1048
mkdir $tmpfolder
cd $data

head -16 $data/${l%.rank}'.res' >> $tmpfolder/resfile
echo 'No.of matches to output .............. 1' >> $tmpfolder/resfile
head -25 $data/${l%.rank}'.res' | tail -8 >> $tmpfolder/resfile

#grep ' H ' $data/$l | head -1 >> $tmpfolder/tophit;
#grep ' M ' $data/$l | head -1 >> $tmpfolder/tophit;
#grep ' A ' $data/$l | head -1 >> $tmpfolder/tophit; 
#grep ' I ' $data/$l | head -1 >> $tmpfolder/tophit;

#if [[ $(head -1 $tmpfolder/tophit) == *"H"* ]]; then
#    quality=H
#elif [[ $(head -1 $tmpfolder/tophit) == *"M"* ]]; then
#    quality=M
#elif [[ $(head -1 $tmpfolder/tophit) == *"A"* ]]; then
#    quality=A
#else
#    quality=I
#fi

awk 'BEGIN{header="n"}{if (substr($1,1,1) ~ "#") header="y"}{if (header ~ "y") print $0} ' $data/$l | grep -v "#" >> $tmpfolder/tophit
pose=$(head -1 $tmpfolder/tophit)							  						##### pose selection
pose=${pose:0:6}
pose="$(echo -e "${pose}" | tr -d '[:space:]')"
if [ ${#pose} == 4 ]; then 
    pose=" "$pose
elif [ ${#pose} == 3 ]; then
    pose="  "$pose
elif [ ${#pose} == 2 ]; then
    pose="   "$pose
elif [ ${#pose} == 1 ]; then
    pose="    "$pose
fi
awk -v p="$pose" '$0 ~ "^" p {print $0}' $data/${l%.rank}'.res' >> $tmpfolder/resfile
pose="$(echo -e "${pose}" | tr -d '[:space:]')"


complex=${l%1-*}
pdb1=$complex'1.pdb'
pdb2=$complex'2.pdb' 

ml singularity/3.5.3
singularity exec --nv /home/g/gabriele/pfs/singularity_img/tf_env.simg python /home/g/gabriele/pfs/utilities/interface_extract.py $data/$pdb1 $data/$pdb2 $tmpfolder/



export LD_LIBRARY_PATH=/home/g/gabriele/pfs/libgsl
export TMDOCKDAT=/home/g/gabriele/pfs/programs/aace_gramm2/data
/home/g/gabriele/pfs/programs/aace_gramm2/aace_gramm -r $tmpfolder/resfile -b true -t AACE18 -d 6.9 -k 5 -o $tmpfolder/${l%1-*}'2' -n 1

mv $tmpfolder/${l%1-*}2_$pose'.pdb' $tmpfolder/${l%1-*}2_$quality
cp $data/$pdb1 $tmpfolder/${pdb1%.pdb}'_'$quality

singularity exec --nv /home/g/gabriele/pfs/singularity_img/tf_env.simg python /home/g/gabriele/pfs/utilities/interface_extract.py $tmpfolder/${pdb1%.pdb}'_'$quality $tmpfolder/${l%1-*}'2_'$quality $tmpfolder/



rm $tmpfolder/${pdb1%.pdb}'_'$quality
rm $tmpfolder/${l%1-*}'2_'$quality
#rm $tmpfolder/resfile
#rm $tmpfolder/tophit

done<$1
