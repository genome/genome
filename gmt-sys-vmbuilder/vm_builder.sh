#!/bin/bash

if [ "$UID" -ne "0" ]
then
  echo "You must be root."
  exit 1
fi

# set up genome repo
wget -O - -q http://repo.gsc.wustl.edu/ubuntu/files/genome-center.asc | sudo apt-key add -
cat <<EOF | sudo tee /etc/apt/sources.list.d/genome.list
deb http://repo.gsc.wustl.edu/ubuntu lucid main
deb http://repo.gsc.wustl.edu/ubuntu lucid-genome-development main
EOF
aptitude update
aptitude install -y git-core
aptitude install -y libur-perl
aptitude install -y libworkflow-perl
aptitude install -y libwebservice-solr-perl libcache-memcached-perl libtest-mockobject-perl bioperl libregexp-common-perl libmime-lite-perl libfile-grep-perl libfile-slurp-perl libinline-perl unzip libdatetime-perl
# refalign
aptitude install -y libfile-copy-recursive-perl

# Somatic

# get oracle-xe-universal_10.2.0.1-1.0_i386.deb
aptitude install -y libaio-dev libaio1

# gmt
aptitude install -y libsort-naturally-perl libanyevent-perl libtest-class-perl libexception-class-perl libmail-sender-perl libguard-perl libcarp-assert-perl libdatetime-perl libfile-ncopy-perl libmail-sendmail-perl libpoe-perl libpoe-component-ikc-perl libtext-csv-perl libmethod-alias-perl libmoosex-types-perl libnamespace-autoclean-perl libmoosex-types-path-class-perl libstatistics-basic-perl libmath-basecalc-perl libgtk2-perl libtest-output-perl libmoosex-strictconstructor-perl liblog-log4perl-perl liblog-dispatch-perl libdbix-dbschema-perl libdbix-class-schema-loader-perl libjson-perl libfile-copy-recursive-perl libgtk2.0-dev libgtk2.0-common cups-client default-jdk pari-gp libgraphics-gnuplotif-perl

aptitude install -y libtext-table-perl libemail-simple-perl libemail-valid-perl unzip weka

aptitude install -y pdl libplack-perl libfile-chdir-perl

# modules to cpan
# FASTAParse
# Class::DBI::Oracle

# modules that need manual install
# do we really need them?
# Math::Pari
# Inline::Java
# 1. Setup JAVA_HOME, then cpan Inline::Java may work

# setup bioperl
cd $HOME
mkdir build
cd $HOME/build
wget ftp://bioperl.org/pub/bioperl/DIST/bioperl-ext-1.5.1.tar.gz
tar bioperl-ext-1.5.1.tar.gz
cd bioperl-ext-1.5.1
perl Makefile.pl

cd $HOME/build
wget http://pari.math.u-bordeaux.fr/pub/pari/unix/pari-2.3.5.tar.gz
tar xvf pari-2.3*
wget http://search.cpan.org/CPAN/authors/id/I/IL/ILYAZ/modules/Math-Pari-2.01080605.tar.gz
tar xvf Math-Pari*
cd Math-Pari*
perl Makefile.pl
make
make install
# could run tests, but it wont make much difference
