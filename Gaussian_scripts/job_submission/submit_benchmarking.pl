#!/usr/bin/perl

use POSIX qw/floor/;

######################
$version=0;
$hours0=3;
$mins0=0;

@l1i=(0); # Basis sets
#@l2i = (0..90); #Distances
@l2i = (1,2); #Distances

@l1a=("3-21G","6-31G","6-311G","cc-pVDZ","cc-pVTZ","cc-pVQZ","cc-pV5Z","cc-pV6Z","CBSB7");
@l2a = (10..100);
$scalar=0.1;
foreach $x (@l2a) { $x = $x * $scalar; }

$l1l='BASS';
$l2l='DIST';

$nMols_per_Task=-1; # -1 to fix the number of cores
$nCores=32;
$nTasks_per_Node=32;

$workDir=`pwd`;
chomp($workDir);

foreach $l1 (@l1i) {

 $l1s=$l1a[$l1];
 $l1d="${l1l}_${l1s}";
 
 # create/change dir
 if (!(-d $l1d)) {
  mkdir($l1d);
 }
 chdir("$workDir/$l1d");

 foreach $l2 (@l2i) {
 
  $l2s=sprintf("%03d",$l2a[$l2]);
  $l2d="${l2l}_${l2s}";
  
  # create/change dir
  if (!(-d $l2d)) {
   mkdir($l2d);
  }
  chdir("$workDir/$l1d/$l2d");

  if ($nMols_per_Task < 0) {
   $nTasks=$nCores;
   $hours_calc=$hours0;
   $mins_calc=$mins0;
  } else {
   $nTasks=$l1a[$l1]/$nMols_per_Task;
   $hours_calc=$hours0;
   $mins_calc=$mins0;
  }

  system "sed -e \'s/##$l1l##/$l1s/g\' -e \'s/##$l2l##/$l2s/g\' $workDir/LiF_CISD.template > LiF_$l1s_$l2s.com";

  # Args[default]: queue[-1] cores cores/node[-1] node-type[-1] mem/core[-1] hours mins name inp out #links exe
  system "submit_Gaussian.pl -1 $nTasks $nTasks_per_Node -1 -1 $hours_calc $mins_calc LiF_$l1s_$l2s LiF_$l1s_$l2s.com LiF_$l1s_$l2s 1 $version";
 
  chdir("$workDir/$l1d");
 
 } #level2

 chdir("$workDir");

} #level1
