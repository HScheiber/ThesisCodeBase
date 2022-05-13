#!/usr/bin/perl

use POSIX qw/floor/;

######################
$cp2k_ver=0;
$hours0=50;
$mins0=0;
#$scaling=1.81;

@l1i=(6);

@l1a=(64,128,256,512,1024,2048,4096);

$l1l='NMOL';

$nMols_per_Task=-1; # -1 to fix the number of cores
$nCores=1024;
$nTasks_per_Node=32;
$out_id="00";
# LOG: 00: 256 cores
######################

$workDir=`pwd`;
chomp($workDir);

foreach $l1 (@l1i) {

 $l1s=sprintf("%05d",$l1a[$l1]);
 $l1d="${l1l}_${l1s}";

 # create/change dir
 if (!(-d $l1d)) {
  mkdir($l1d);
 }
 chdir("$workDir/$l1d");

 if ($nMols_per_Task < 0) {
  $nTasks=$nCores;
  $hours_calc=$hours0;
  $mins_calc=$mins0;
  #$mins_calc=floor(((2**$l1)**$scaling)*($hours0*60+$mins0));
 } else {
  $nTasks=$l1a[$l1]/$nMols_per_Task;
  $hours_calc=$hours0;
  $mins_calc=$mins0;
 }

 $ids="$l1d";
 $sTasks=sprintf("%05d",$nTasks);

 $path = "$workDir";
 $path =~ s/\//\\\//g;

 system "sed -e \'s/##$l1l##/$path\\\/templates-liq\\\/$l1s/g\' $workDir/input.template > input.$out_id.inp";
  system "echo Submitting $l1s Molecules at $nCores cores for $mins_calc mins";

 # Args[default]: queue[-1] cores cores/node[-1] node-type[-1] mem/core[-1] hours mins name inp out #links exe
 system "submit_cp2k.pl -1 $nTasks $nTasks_per_Node -1 -1 $hours_calc $mins_calc StrScOT-$l1s-$nCores input.$out_id.inp output.$out_id.$sTasks.out 1 $cp2k_ver";
  
 chdir("$workDir");

} #level1

