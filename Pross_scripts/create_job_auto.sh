#!/bin/sh

#produces one pair of files, commandline and jobname. Commandline is used to execute jobname through lsf on the fleishman queue.
#USAGE: for i in *.pdb; do fleish_sub.sh -s ${i} @flags; done

dir=$1 
date=`date +%s%N | awk '{print substr($0,9,length($0)-11)}'`
jobname=${dir}/job/job.${date}
commandname=${dir}/command
if [ ! -d ${dir}/job ];then
	mkdir ${dir}/job
fi
if [ ! -d ${dir}/out ];then
	mkdir ${dir}/out
fi
if [ ! -d ${dir}/err ];then
	mkdir ${dir}/err
fi
if [ ! -d ${dir}/checkpoint ];then
	mkdir ${dir}/checkpoint
fi
if [ ! -d ${dir}/pdbs ];then
	mkdir ${dir}/pdbs
fi
if [ ! -d ${dir}/scores ];then
	mkdir ${dir}/scores
fi

if [ ! -d ${dir}/resfiles ];then
        mkdir ${dir}/resfiles
fi
out=${dir}/out/out.${date}
err=${dir}/err/err.${date}
check=${dir}/checkpoint/

echo "#!/bin/bash" > $jobname
#echo ". /usr/share/lsf/conf/profile.lsf" >> $jobname //Jerry Mershel told me to take this line off (email from 18/9)

echo "cd $dir">> $jobname
printf "/home/labs/fleishman/adig/Rosetta/main/source/bin/rosetta_scripts.static.linuxgccrelease" >> $jobname

until [ -z $2 ]; do
	printf " $2 "
	shift
done >> $jobname

echo "" >> $jobname

echo "bsub -u /dev/null  -R rusage[mem=2048] -L /bin/bash -G fleishman-wx-grp-lsf -q new-all.q -o $out -e $err /apps/RH6U4/blcr/0.8.5/bin/cr_run $jobname" >> $commandname
#echo "bsub -u /dev/null  -R rusage[mem=2048] -L /bin/bash -G fleishman-wx-grp-lsf -q fleishman -o $out -e $err /apps/RH6U4/blcr/0.8.5/bin/cr_run $jobname" >> $commandname
#echo "bsub -u /dev/null -L /bin/bash fleishman-grp-lsf -q new-all.q -o $out -e $err -k $check $jobname" >> $commandname 
chmod +x $jobname

