############## creat and submit multiple jobs
#!/bin/sh
N_job=10
workdir=`pwd`

for((k=1;k<=$N_job;k++))
do
    cd $workdir
	cd ..

    if [ ! -d data_$k ]; then
        echo "No data_${k}"
    else	
	    #rm -r data_${k}/output_hadron 	
        #cp -r ./output_hadron data_${k}
		
	    rm -r data_${k}/output_hadron/positive/frag 	
        cp -r ./output_hadron/positive/frag data_${k}/output_hadron/positive		
		
        cd data_${k}/output_hadron/positive/frag
        ./exec_frag.sh 1>log 2>err &
        #cd data_${k}/output_hadron/negative/frag
        #./exec_frag.sh 1>log 2>err &
    fi
done
