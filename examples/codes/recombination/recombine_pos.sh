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
		
	    rm -r data_${k}/output_hadron/positive 	
        cp -r ./output_hadron/positive data_${k}/output_hadron/		
		
        cd data_${k}/output_hadron/positive
        ./exec.sh 1>log 2>err &
        #cd data_${k}/output_hadron/positive
        #./exec.sh 1>log 2>err &
    fi
done
