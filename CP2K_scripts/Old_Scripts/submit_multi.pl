#!/usr/bin/perl

use POSIX qw/floor/;

######################
$cp2k_ver=0;
$hours0=8;
$mins0=0;
$scaling=0.0;

@l1i=(1);
@l2i=(3);
@l3i=(0,1,2,3,4,5,6,7,8);

@l1a=(1.2,1.6,2.0);
@l2a=(1E-5,1E-4,1E-3,1E-2);
@l3a=(5E-5,1E-5,5E-6,1E-6,0,-1E-6,-5E-6,-1E-5,-5E-5);

$l1l='XARC';
$l2l='EPEE';
$l3l='SHGA';

$nMols_per_Task=-1;
$nCores=32;
$nTasks_per_Node=16;
$out_id="01";
# LOG: 00: 512 cores
#      01: 1024 cores
######################

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
 
  $l2s=$l2a[$l2];
  $l2d="${l2l}_${l2s}";
 
  # create/change dir
  if (!(-d $l2d)) {
   mkdir($l2d);
  }
  chdir("$workDir/$l1d/$l2d");
 
  foreach $l3 (@l3i) {
 
   $l3s=$l3a[$l3];
   $l3d="${l3l}_${l3s}";
  
   # create/change dir
   if (!(-d $l3d)) {
    mkdir($l3d);
   }
   chdir("$workDir/$l1d/$l2d/$l3d");

   if ($nMols_per_Task < 0) {
    $nTasks=$nCores;
    $hours_calc=0;
    $mins_calc=floor(((2**$l1)**$scaling)*($hours0*60+$mins0));
   } else {
    $nTasks=$l1a[$l1]/$nMols_per_Task;
    $hours_calc=$hours0;
    $mins_calc=$mins0;
   }

   $ids="$l1d-$l2d-$l3d";
   $sTasks=sprintf("%05d",$nTasks);

   system "sed -e \'s/##$l1l##/$l1s/g\' -e \'s/##$l2l##/$l2s/g\' -e \'s/##$l3l##/$l3s/g\' $workDir/input.template > input.$out_id.inp";

   # Args[default]: queue[-1] cores cores/node[-1] node-type[-1] mem/core[-1] hours mins name inp out #links exe
   system "submit_cp2k_new.pl -1 $nTasks $nTasks_per_Node sandybridge -1 $hours_calc $mins_calc b${l1}${l2}${l3}-$sTasks input.$out_id.inp output.$out_id.$sTasks.out 4 $cp2k_ver";
  
   chdir("$workDir/$l1d/$l2d");
 
  } #level3
 
  chdir("$workDir/$l1d");
 
 } #level2

 chdir("$workDir");

} #level1
