#!/bin/bash
SCRIPT_PATH=`readlink -f $0`
DIR_PATH=`dirname $SCRIPT_PATH`
# DIR_PATH is the location of the vmbuilder directory
GENOME_DIR=`dirname $DIR_PATH`

if [ "$UID" -eq "0" ]
then
  echo "This must be run as a regular user"
  exit 1
fi

# source lsf file if it exists
cat <<EOF | tee -a $HOME/.bashrc
if [ -f /usr/local/lsf/conf/profile.lsf ]; then
    . /usr/local/lsf/conf/profile.lsf
fi

EOF

echo "export PERL5LIB=$GENOME_DIR/lib/perl:\$PERL5LIB" >> $HOME/.bashrc
echo "export PATH=$GENOME_DIR/bin:\$PATH" >> $HOME/.bashrc
echo "export TNS_ADMIN=$DIR_PATH" >> $HOME/.bashrc
echo "export ORACLE_HOME=/opt/oracle-instantclient" >> $HOME/.bashrc

