#!/usr/bin/perl

use POSIX qw(ceil floor);

$args=13;
if ($#ARGV!=$args) {
 print "Incorrect number of inputs.\n";
 die "Args[default]: queue[-1] cores cores/node[-1] node-type[-1] mem/core[-1] hours mins name inp out \n";
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
 $submit=$ARGV[13];
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
 $selected_exe="runcry17";
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

for $i (1..$numlinks) {

  if ($i == 1 and $restart_switch == 1) {
  $CompTag="RUN.COMPLETE";
  ## Remove any old COMPLETE.RUN files if forcing a job restart
    if(-e $CompTag) {
      unlink ($CompTag);
    }
  }


  $strI=sprintf("%03d", $i);
  $linkmainname="$mainname-$strI";
  print "CREATING CRYSTAL17 BATCH SCRIPT FOR LINK: $i\n";
  open(SUBM,">$linkmainname.subm");

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
#PBS -N $linkmainname
#PBS -e $linkmainname.stde
#PBS -o $linkmainname.stdo

VERBATIM

  }
  elsif ($scheduler==1) { #SLURM

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
#SBATCH --error=_$mainname.stde
#SBATCH --export=ALL

VERBATIM

  }

  print SUBM <<"VERBATIM";
# set EXE environment
PATH=\$PATH
export PATH
EXE=$selected_exe

VERBATIM
    
  if ($scheduler==2) { #PBS
    print SUBM <<"VERBATIM";
export CRY17_SCRDIR=/global/scratch/\$USER/\$PBS_JOBID

# Preparing and moving inputfiles to tmp:
submitdir=\$PBS_O_WORKDIR
cd \$submitdir
VERBATIM
  }
  elsif ($scheduler==1) { #SLURM
    print SUBM <<"VERBATIM";
export CRY17_SCRDIR=/scratch/\$USER/\$SLURM_JOB_ID

# Preparing and moving inputfiles to tmp:
submitdir=\$SLURM_SUBMIT_DIR
cd \$submitdir
VERBATIM
 }

if ($i != 1) {
  if ($scheduler==2) { #PBS
    print SUBM <<"VERBATIM";
module load python/3.5.0
pyRS=`which Restart_CRYSTAL.py`
python3 \$pyRS $mainname 1 > /dev/null
value=\$(<Restart.Response)
rm Restart.Response

if [[ "\$value" == "0" ]];
then
    exit 0
fi

VERBATIM
  }
  
  elsif ($scheduler==1) { #SLURM
    print SUBM <<"VERBATIM"; 
module load python/3.7.0
pyRS=`which Restart_CRYSTAL.py`
python3 \$pyRS $mainname 1 > /dev/null
value=\$(<Restart.Response)
rm Restart.Response

if [[ "\$value" == "0" ]];
then
    exit 0
fi

VERBATIM
  }
  
}
elsif ($i == 1 and $restart_switch == 1){
  if ($scheduler==2) { #PBS
    print SUBM <<"VERBATIM";
module load python/3.5.0
pyRS=`which Restart_CRYSTAL.py`
python3 \$pyRS $mainname $restart_WF > /dev/null
value=\$(<Restart.Response)
rm Restart.Response

if [[ "\$value" == "0" ]];
then
    exit 0
fi

VERBATIM
  }
  
  elsif ($scheduler==1) { #SLURM
    print SUBM <<"VERBATIM"; 
module load python/3.7.0
pyRS=`which Restart_CRYSTAL.py`
python3 \$pyRS $mainname $restart_WF > /dev/null
value=\$(<Restart.Response)
rm Restart.Response

if [[ "\$value" == "0" ]];
then
    exit 0
fi

VERBATIM
  }
}
 
 print SUBM <<"VERBATIM";

######################################
# Section for running the program:

# Running the program:
\$EXE $name_in $name_in

echo "Job finished at"
date
pyCC=`which Check_Complete.py`
python3 \$pyCC $mainname > /dev/null
################### Job Ended ###################
exit 0

VERBATIM

close(SUBM);
 if ( $submit==1) {
  if ($i==1) {
    print "SUBMITTING LINK: 1.\n\n";
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
 } else {
  print "LINK: $i BATCH SCRIPT GENERATED BUT NOT SUBMITTED.\n\n";
 }
}


