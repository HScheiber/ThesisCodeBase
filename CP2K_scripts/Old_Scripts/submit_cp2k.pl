#!/usr/bin/perl

use POSIX qw(ceil floor);

# globals
$scheduler=1; # 2 - PBS, 1 - SLURM
#$account="mfj-281-ab"; # charging account
$account="def-rzk"; # charging account
@exe = (
 "cp2k.popt",
 "/project/6004141/group/cp2k/cp2k-mcgill/cp2k/exe/Graham-Linux-x86-64-intel-intelmpi/cp2k.popt"
);

$args=11;
if ($#ARGV!=$args) {
 print "Versions of executable. The first option - CP2K in your path:\n";
 for ($ii=0; $ii<=$#exe; $ii++) {
  printf "%3d - %-100s\n", $ii, $exe[$ii];
 }
 die "Args[default]: queue[-1] cores cores/node[-1] node-type[-1] mem/core[-1] hours mins name inp out #links exe\n";
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
 if ($ARGV[11]>$#exe || $ARGV[11]<0) {die "Check exe number";}
 $selected_exe=$exe[$ARGV[11]];
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
#PBS -o _$linkmainname.stdo

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
      print SUBM "#SBATCH --mem=MaxMemPerNode\n";
	}
 
    print SUBM <<"VERBATIM";
#SBATCH --account=$account 
#SBATCH --job-name=$linkmainname
#SBATCH --error=_$linkmainname.stde
#SBATCH --output=_$linkmainname.stdo

VERBATIM

  }

  print SUBM <<"VERBATIM";
# set EXE environment
PATH=\$PATH:\$HOME/ThesisCodeBase
export PATH
VERBATIM
    
  print SUBM "EXE=$selected_exe\n";
  if ($scheduler==2) {
    print SUBM "workd=\${PBS_O_WORKDIR}\n";
  }
  elsif ($scheduler==1) {
    print SUBM "workd=\${SLURM_SUBMIT_DIR}\n";
  }
  print SUBM "cd \$workd\n";

  if ($i != 1) {
    $name_in=$mainname."-1.restart";
  }
  #print SUBM "mpirun -n $cores \$EXE $name_in > \$workd/$name_out-$strI\n";
  print SUBM "srun \$EXE $name_in > \$workd/$name_out-$strI\n";

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


