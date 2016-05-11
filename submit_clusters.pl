#!/usr/bin/perl

use POSIX qw(ceil floor);

# globals
$scheduler=2; # 2 - PBS, 1 - SLURM
$account="mfj-281-aa"; # charging account

$args=9;
if ($#ARGV!=$args) {
 die "Args: queue-name[-1] cores cores-per-node[-1] node-type[-1] mem-per-core[-1] hours mins name inp out\n";
}
else {
 $queuename=$ARGV[0];
 $cores=$ARGV[1];
 $pernode=$ARGV[2];
 $nodetype=$ARGV[3];
 $mempercore=$ARGV[4];
 $hours=$ARGV[5];
 $mins=$ARGV[6];
 $mainname=$ARGV[7];
 $name_in=$ARGV[8];
 $name_out=$ARGV[9];
}

# default values
# number of nodes
if ($pernode>0) {
 $nodes=ceil($cores/$pernode);
}
else {
 $node=-1;
}
# type of nodes: used only if $pernode is set
if ($pernode>0) {
 # check if the value is allowed
 if ($nodetype =~ /westmere/) {
  $nodetype = "westmere";
 }
 elsif ($nodetype =~ /sandybridge/) {
  $nodetype = "sandybridge";
 }
 elsif ($nodetype =~ /-1/) {
  $nodetype = "-1";
 }
 else {
  print "Illegal value node type\n";
  exit 1;
 }
}
else {
 $nodetype="-1";
}

print "CREATING BATCH SCRIPT\n";
open(SUBM,">$mainname.subm");

if ($scheduler==1) { # SLURM
 
 print SUBM <<"VERBATIM";
#!/bin/csh
#SBATCH --time=$hours:$mins:00
#SBATCH --ntasks=$cores
VERBATIM
 
 if ($pernode>0) {
  print SUBM "#SBATCH --ntasks-per-node=$pernode\n";
  #print SUBM "#SBATCH --mem-per-cpu=$mempercore\n";
 }
 print SUBM "#SBATCH --account=$account\n";

 print SUBM <<"VERBATIM";
#SBATCH --job-name=$mainname
#SBATCH --error=_$mainname.stde
#SBATCH --output=_$mainname.stdo

# set EXE environment
set path = (\$path ~/bin)
#set EXE = \$HOME/bin/cp2k.popt
set EXE = cp2k.popt

# WORKING directory
set workd = \${SLURM_SUBMIT_DIR}

# some scratch related variables
set scrroot = \$SCRATCH/$mainname.\${SLURM_JOB_ID}

cd \$workd

# set the SCRATCH directory
mkdir -p \$scrroot

# bring input files to the scratch dir
cp $name_in \$scrroot
cp * \$scrroot
if ( -e coords.i ) then
 cp coords.i \$scrroot
endif

cd \$scrroot
VERBATIM

 if ($pernode>0) {
  print SUBM "aprun -n $cores -N $pernode \$EXE $name_in > \$workd/$name_out\n";
 }
 else {
  print SUBM "aprun -n $cores \$EXE $name_in > \$workd/$name_out\n";
 }
 print SUBM <<"VERBATIM";

# get output files
tar -zcvf $name_out.tar.gz *
mv $name_out.tar.gz \$workd/

cd \$workd

/bin/rm -rf \$scrroot

VERBATIM

} elsif ($scheduler==2) { # PBS

 print SUBM <<"VERBATIM";
#!/bin/bash
#PBS -l walltime=$hours:$mins:00
VERBATIM

 if ($pernode>0) {
  if ($nodetype =~ "-1") {
   print SUBM "#PBS -l nodes=$nodes:ppn=$pernode\n";
  }
  else {
   print SUBM "#PBS -l nodes=$nodes:ppn=$pernode:$nodetype\n";
  }
 }
 else {
  print SUBM "#PBS -l procs=$cores\n";
 }

 if ($mempercore>0) {
  print SUBM "#PBS -l pmem=$mempercore\n";
 }
 
 print SUBM <<"VERBATIM";
#PBS -A $account 
#PBS -V
#PBS -N $mainname
#PBS -e _$mainname.stde
#PBS -o _$mainname.stdo

VERBATIM

  print SUBM <<"VERBATIM";
# set EXE environment
PATH=\$PATH:\$HOME/bin
export PATH
EXE=\$HOME/bin/cp2k.popt

# WORKING directory
workd=\${PBS_O_WORKDIR}

cd \$workd

tar -zxf packed.tar.gz

VERBATIM

 print SUBM "mpirun -n $cores \$EXE $name_in > \$workd/$name_out\n";
 print SUBM <<"VERBATIM";

rm packed.tar.gz
tar --remove-files -zcf packed.tar.gz *

VERBATIM

} else {
 print "Unknown scheduler type";
 exit 1;
}

close(SUBM);

if ($scheduler==1) {
 system("sbatch $mainname.subm");
} 
elsif ($scheduler==2) {
 if ($queuename==-1) {
  system("qsub $mainname.subm");
 }
 else {
  system("qsub -q $queuename $mainname.subm");
 }
}

