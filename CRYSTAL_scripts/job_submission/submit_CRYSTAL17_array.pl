#!/usr/bin/perl

use POSIX qw(ceil floor);

$args=14;
if ($#ARGV!=$args) {
 print "Incorrect number of inputs.\n";
 die "Args[default]: queue[-1] cores cores/node[-1] node-type[-1] mem/core[-1] hours mins name inp out a-parameters \n";
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
 $numlinks=$ARGV[10];
 $restart_switch=$ARGV[11];
 $restart_WF=$ARGV[12];
 $aparmeters=$ARGV[13];
 $submit=$ARGV[14];
}

$server=`uname -n | cut -c 1-3`;
chomp($server);

if ($server eq 'gra'){
 $account="def-patey"; # charging account for graham
 $scheduler=1; # SLURM scheduler
 $selected_exe="runcry17";
} elsif ($server eq 'sea') {
 $account="rrg-patey-ad"; # charging account for orcinus
 $scheduler=2; # PBS scheduler
 $selected_exe="/home/scheiber/CRYSTAL17/utils17/runcry17";
} elsif ($server eq 'ced') {
 $account="rrg-patey-ad"; # charging account for cedar
 $scheduler=1; # SLURM scheduler
 $selected_exe="runcry17";
} else {
 die "Unknown machine, stopped";
}

if ($scheduler < 1 && $scheduler > 2) {
  print "Unknown scheduler type";
  exit 1;
}

# default values
# number of nodes
if ($pernode>0) {
 $nodes=ceil($cores/$pernode);
}
else {
 $nodes=-1;
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

print "CREATING CRYSTAL17 ARRAY BATCH SCRIPT\n";
open(SUBM,">$mainname.subm");

  if ($scheduler==2) { # PBS
    print SUBM <<"VERBATIM";
#!/bin/bash
#PBS -S /bin/bash
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
    elsif ($mempercore==0) {
      print SUBM "#PBS -l pmem=2000mb\n";
    }
 
    print SUBM <<"VERBATIM";
#PBS -V
#PBS -N $mainname
#PBS -t $aparmeters
#PBS -e $mainname.stde
#PBS -o $mainname.stdo

VERBATIM

  }
  elsif ($scheduler==1) { # SLURM

    print SUBM <<"VERBATIM";
#!/bin/bash
#SBATCH --time=$hours:$mins:00
VERBATIM

    if ($pernode>0) {
      print SUBM "#SBATCH --nodes=$nodes\n";
      print SUBM "#SBATCH --tasks-per-node=$pernode\n";
    }
    else {
      print SUBM "#SBATCH --ntasks=$cores\n";
    }

    if ($mempercore>0) {
      print SUBM "#SBATCH --mem=$mempercore\n";
    }
	elsif ($mempercore==0) {
      print SUBM "#SBATCH --mem=0\n";
	}
 
    print SUBM <<"VERBATIM";
#SBATCH --account=$account 
#SBATCH --job-name=$mainname
#SBATCH --array=$aparmeters
#SBATCH --export=ALL
#SBATCH --output=$mainname.%a.out

VERBATIM

  }

  print SUBM <<"VERBATIM";
# set EXE environment
PATH=\$PATH
export PATH
EXE=$selected_exe

VERBATIM
  
  if ($scheduler==2) { # PBS
    print SUBM <<"VERBATIM";
id=\${PBS_ARRAYID}
anum=\$(echo \"scale=10; \$id / 100;\" | bc)
aparam=\$(printf '%05.2f' \$anum)
idDir=\$(printf 'APAR_%05.2f' \$aparam)

# Preparing and moving inputfiles to tmp:
submitdir=\$PBS_O_WORKDIR
export CRY17_INP=\$submitdir/\$idDir
cd \$submitdir/\$idDir

######################################
# Section for running the program and cleaning up:

# Running the program: 
time \$EXE $name_in-\$aparam

# Cleaning up and moving files back to home/submitdir:
# Make sure to move all essential files specific for the given job/software.

# Clean up.
rm $name_in-\$aparam.f98
rm $name_in-\$aparam.f9
cd \$submitdir
VERBATIM
  }
  elsif ($scheduler==1) { #SLURM
 print SUBM <<"VERBATIM";
id=\${SLURM_ARRAY_TASK_ID}
anum=\$(echo \"scale=10; \$id / 100;\" | bc)
aparam=\$(printf '%05.2f' \$anum)
idDir=\$(printf 'APAR_%05.2f' \$aparam)

export CRY17_SCRDIR=/scratch/\$USER/\$SLURM_ARRAY_JOB_ID-\$idDir
tempdir=\$CRY17_SCRDIR
mkdir -p \$tempdir

# Preparing and moving inputfiles to tmp:
submitdir=\$SLURM_SUBMIT_DIR
cp -r \$submitdir/\$idDir/$name_in-\$aparam.d12 \$tempdir/$name_in-\$aparam.d12
cd \$tempdir

scontrol update jobid=\$SLURM_ARRAY_JOB_ID\\_\$SLURM_ARRAY_TASK_ID jobname=$mainname\_\$aparam

######################################
# Section for running the program and cleaning up:

# Running the program: 
time \$EXE $name_in-\$aparam

# Cleaning up and moving files back to home/submitdir:
# Make sure to move all essential files specific for the given job/software.

cp $name_in-\$aparam.out \$submitdir/\$idDir/$name_out-\$aparam.out

# Clean up.
cd \$submitdir
rm -r \$tempdir/*
rmdir \$tempdir
VERBATIM
  }
  
 print SUBM <<"VERBATIM";

echo "Job finished at"
date
################### Job Ended ###################
exit 0

VERBATIM

close(SUBM);

if ($submit==1) {
 print "SUBMITTING THE CRYSTAL17 ARRAY JOB.\n\n";
 if ($scheduler==2) {
   $FULLJOBID=`qsub $mainname.subm`;
 }
 elsif ($scheduler==1) {
   $FULLJOBID=`sbatch --export=all $mainname.subm`;
 }
} else {
print "CRYSTAL17 JOB SCRIPT GENERATED BUT NOT SUBMITTED.\n\n";
}




