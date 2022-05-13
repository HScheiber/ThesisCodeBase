#!/usr/bin/perl

use POSIX qw(ceil floor);

$args=12;
if ($#ARGV!=$args) {
 print "Incorrect number of inputs.\n";
 die "Args[default]: queue[-1] cores cores/node[-1] node-type[-1] mem/core[-1] hours mins name inp out #links a-parameters properties\n";
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
 $aparmeters=$ARGV[11];
 $properties=$ARGV[12]; # Must be one of: -1 -> Normal; 0 -> DSDPBE86; 1 -> B2GPPLYP
}

$server=`uname -n | cut -c 1-3`;
chomp($server);

if ($server eq 'gra'){
 $account="def-patey"; # charging account for graham
 $scheduler=1; # SLURM scheduler
 if ($properties<0){
  $module="gaussian/g16.b01";
  $selected_exe="g16";
 } elsif ($properties==0 || $properties==1){
  $module="gaussian/g09.e01";
  $selected_exe="g09";
 } else {
  die "Unknown properties input, stopped";
 }
}
elsif ($server eq 'sea') {
 $account="rrg-patey-ad"; # charging account for orcinus
 $scheduler=2; # PBS scheduler
 $module="gaussian/g09d01";
 $selected_exe="g09";
}
elsif ($server eq 'ced') {
 $account="rrg-patey-ad"; # charging account for cedar
 $scheduler=1; # SLURM scheduler
 if ($properties<0){
  $module="gaussian/g16.b01";
  $selected_exe="g16";
 } elsif ($properties==0 || $properties==1){
  $module="gaussian/g09.e01";
  $selected_exe="g09";
 } else {
  die "Unknown properties input, stopped";
 }
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

for $i (1..$numlinks) {
  $strI=sprintf("%03d", $i);
  $linkmainname="$mainname-$strI";
  print "CREATING BATCH SCRIPT FOR LINK: $i\n";
  open(SUBM,">$linkmainname.subm");

  if ($scheduler==2) { # PBS
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
#PBS -N $linkmainname
#PBS -e _$linkmainname.stde

VERBATIM

  }
  elsif ($scheduler==1) {

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

    #if ($nodetype =~ "-1") {
    #  print SUBM "#SBATCH $nodetype\n";
    #}

    if ($mempercore>0) {
      print SUBM "#SBATCH --mem=$mempercore\n";
    }
	else {
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
module load $module
PATH=\$PATH
export PATH
VERBATIM

if ($properties==0){
 print SUBM <<"VERBATIM";
export GAUSS_DFTD3_S6=480000
export GAUSS_DFTD3_SR6=0
export GAUSS_DFTD3_S8=0
export GAUSS_DFTD3_ABJ1=0
export GAUSS_DFTD3_ABJ2=5600000
VERBATIM
}
elsif ($properties==1){
 print SUBM <<"VERBATIM";
export GAUSS_DFTD3_S6=560000
export GAUSS_DFTD3_SR6=0
export GAUSS_DFTD3_S8=259700
export GAUSS_DFTD3_ABJ1=0
export GAUSS_DFTD3_ABJ2=6333200
VERBATIM
}
    
  print SUBM "EXE=$selected_exe\n";
  if ($scheduler==2) {
    print SUBM "workd=\${PBS_O_WORKDIR}\n";
  }
  elsif ($scheduler==1) {
 print SUBM <<"VERBATIM";
id=\${SLURM_ARRAY_TASK_ID}
anum=\$(echo \"scale=10; \$id / 100;\" | bc)
aparam=\$(printf '%05.2f' \$anum)
idDir=\$(printf 'APAR_%05.2f' \$aparam)
VERBATIM
  }
  
 print SUBM <<"VERBATIM";

LINDA=\${SLURM_JOB_NODELIST}
corename=\$(echo \$LINDA | cut -c 1-3)
corelist=\$(echo \$LINDA | sed -e "s/,/,\$corename/g" -e "s/\\[//g" -e "s/\\]//g")

export GAUSS_SCRDIR=/scratch/\$USER/\$SLURM_ARRAY_JOB_ID/\$id
tempdir=\$GAUSS_SCRDIR
mkdir -p \$tempdir

# Preparing and moving inputfiles to tmp:
submitdir=\$SLURM_SUBMIT_DIR
sed -e "s/##NODES##/\$corelist/g" \$submitdir/\$idDir/$name_in-\$aparam.com > \$tempdir/$name_in-\$aparam.com
cd \$tempdir

scontrol update jobid=\$SLURM_ARRAY_JOB_ID\\_\$SLURM_ARRAY_TASK_ID jobname=$mainname\_\$aparam

######################################
# Section for running the program and cleaning up:

# Running the program:
VERBATIM
  
  if ($i != 1) {
    $name_in=$mainname."-1.restart";
  }

 print SUBM <<"VERBATIM";
 
time \$EXE < $name_in-\$aparam.com >& $name_out-\$aparam.log

# Cleaning up and moving files back to home/submitdir:
# Make sure to move all essential files specific for the given job/software.

cp $name_out-\$aparam.log \$submitdir/\$idDir/$name_out-\$aparam.log

# Clean up.
cd \$submitdir
rm  \$tempdir/*
rmdir \$tempdir

echo "Job finished at"
date
################### Job Ended ###################
exit 0

VERBATIM

  close(SUBM);
  
  if ($i==1) {
    print "SUBMITTING THE FIRST LINK $i.\n\n";
    if ($scheduler==2) {
      $FULLJOBID=`qsub $linkmainname.subm`;
    }
    elsif ($scheduler==1) {
      $FULLJOBID=`sbatch --export=all $linkmainname.subm`;
    }
  }
  else {
    print "SUBMITTING LINK: $i.\n\n\n";
    if ($scheduler==2) {
      @SPLITJOBID= split /\./,$FULLJOBID;
      $JOBID=$SPLITJOBID[0];
      $FULLJOBID=`qsub -W depend=afterany:$JOBID $linkmainname.subm`;
    }
    elsif ($scheduler==1) {
      @SPLITJOBID= split /\s+/,$FULLJOBID;
      $JOBID=$SPLITJOBID[3];
      print "$JOBID\n";
      $FULLJOBID=`sbatch --export=all --depend=afterany:$JOBID $linkmainname.subm`;
    }
  }

}



