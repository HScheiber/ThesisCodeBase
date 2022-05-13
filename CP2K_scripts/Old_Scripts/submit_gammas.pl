#!/usr/bin/perl

use POSIX qw/floor/;

######################
$cp2k_ver=0;
$hours0=10;
$mins0=0;
$scaling=0.0;

@l1i=(0,1,2,3,4);
@l1a=(1E-5,5E-5,5E-4,1E-3,5E-3);

$l1l='SHGA';
$l2l='NAME';

$nMols_per_Task=-1;
$nCores=64;
$nTasks_per_Node=16;
$out_id="01";
# LOG: 00: OT+Langevin
# LOG: 01: ALMO_APP+Langevin
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
 
 if ($nMols_per_Task < 0) {
  $nTasks=$nCores;
  $hours_calc=0;
  $mins_calc=floor(((2**$l1)**$scaling)*($hours0*60+$mins0));
 } else {
  $nTasks=$l1a[$l1]/$nMols_per_Task;
  $hours_calc=$hours0;
  $mins_calc=$mins0;
 }

 $sTasks=sprintf("%05d",$nTasks);

 system "sed -e \'s/##$l1l##/$l1s/g\' -e \'s/##$l2l##/b${l1}-$sTasks/g\' $workDir/input.template > input.$out_id.inp";

 # Args[default]: queue[-1] cores cores/node[-1] node-type[-1] mem/core[-1] hours mins name inp out #links exe
 system "submit_cp2k_new.pl -1 $nTasks $nTasks_per_Node sandybridge -1 $hours_calc $mins_calc b${l1}-$sTasks input.$out_id.inp output.$out_id.$sTasks.out 1 $cp2k_ver";

 chdir("$workDir");

} #level1
