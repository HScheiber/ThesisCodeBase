#!/usr/bin/perl

use POSIX qw/floor/;

######################
$cp2k_ver=0;
$hours0=8;
$mins0=0;
$scaling=0;

@l1i=(11);
@l2i=(3);

@l1a=(64,128,256,512,1024,2048,4096,6144,8192,12288,16384,32768);
@l2a=(1,2,4,8,16);
$l3s=10;

$l1l='NMOL';
$l2l='MCOR';
$l3l='MTAR';

#$nMols_per_Task=-1; # -1 to fix the number of cores
#$nCores=256;
#$nTasks_per_Node=16;
$out_id="00";
# LOG: 00: 16 molecules per core
# LOG: 01: 8 molecules per core
# LOG: 02: 4 molecules per core
# LOG: 03: 2 molecules per core
# LOG: 04: 1 molecules per core
######################

$workDir=`pwd`;
chomp($workDir);

foreach $l1 (@l1i) {

 $l1s=sprintf("%05d",$l1a[$l1]);
 $l1d="${l1l}_${l1s}";
 
 # Set matrix iterate target factor higher for 32k system
 # if ($l1a[$l1] < 16500) {
  # $l3s=50;
 # } else {
  # $l3s=1000;
 # }

 
 # create/change dir
 if (!(-d $l1d)) {
  mkdir($l1d);
 }
 chdir("$workDir/$l1d");

 foreach $l2 (@l2i) {
 
  $l2s=$l2a[$l2];
  $nMols_per_Task = $l2a[$l2];
  $l2d="${l2l}_${l2s}";
 
  # create/change dir
  if (!(-d $l2d)) {
   mkdir($l2d);
  }
  chdir("$workDir/$l1d/$l2d");
 
   if ($nMols_per_Task < 0) {
    $nTasks=$nCores;
    $hours_calc=0;
    $mins_calc=floor(((2**$l1)**$scaling)*($hours0*60+$mins0));
   } else {
    $nTasks=$l1a[$l1]/$nMols_per_Task;
    #$hours_calc=($hours0*$l2s);
    #$mins_calc=($mins0*$l2s);
	$hours_calc=$hours0;
	$mins_calc=$mins0;
   }

   if ($nTasks < 32) {
    $nTasks_per_Node=$nTasks;
   } else {
    $nTasks_per_Node = 32;
   }
   
   $ids="$l1d-$l2d";
   $sTasks=sprintf("%05d",$nTasks);

   $path = "$workDir";
   $path =~ s/\//\\\//g;

   system "sed -e \'s/##$l1l##/$path\\\/templates-liq\\\/$l1s/g\' -e \'s/##$l3l##/$l3s/g\' $workDir/input.template > input.$out_id.inp";
   system "echo Submitting $l1s Molecules at $l2s Molecules per core for wall time $hours_calc:$mins_calc";

   # Args[default]: queue[-1] cores cores/node[-1] node-type[-1] mem/core[-1] hours mins name inp out #links exe
   system "submit_cp2k.pl -1 $nTasks $nTasks_per_Node -1 -1 $hours_calc $mins_calc $l2s-${l1}${l2}-$l1s input.$out_id.inp output.$out_id.$sTasks.out 1 $cp2k_ver";
  
  chdir("$workDir/$l1d");
 
 } #level2

 chdir("$workDir");

} #level1

