#!/gsc/bin/sh
ID=$1;
DIR=$2;
cd ${DIR};
    g=`cut -f3,5 $ID.cfg | grep normal | cut -f2 | perl -ane 'chomp; @x = split ":"; print "$x[1]\n"' | sort -u`;
    echo $g;
    g=`echo $g | perl -ane 'chomp; $_=~s/ /,/g; print $_'`
    if [ -e $ID.novo.ctx ]; then 
	echo ./AssembleSV_chr_new_LIB.sh $ID $g .novo ${DIR}; 
    fi; 
