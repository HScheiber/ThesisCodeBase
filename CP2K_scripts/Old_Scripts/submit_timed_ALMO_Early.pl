#!/usr/bin/perl

use POSIX qw/floor/;
use warnings;

######################
$cp2k_ver=0;
$hours0=3;
$mins0=0;
$scaling=1.0;

@l1i=(0,1,2);
@l2i=(0,1,2,3,4,5,6);


@l1a=(1E-2,1E-3,1E-4);
@l2a=(64,128,256,512,1024,2048,4096,8192);

$l1l='EPEE';
$l2l='NMOL';

$nMols_per_Task=-1; # -1 to fix the number of cores
$nCores=1024;
$nTasks_per_Node=32;
$out_id="00";
# LOG: 00: 1048 cores
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

  $l2s=sprintf("%05d",$l2a[$l2]);
  $l2d="${l2l}_${l2s}";

  # create/change dir
  if (!(-d $l2d)) {
  mkdir($l2d);
  }
  chdir("$workDir/$l1d/$l2d");
 
  if ($nMols_per_Task < 0) {
  $nTasks=$nCores;
  $hours_calc=$hours0;
  #$mins_calc=floor(((2**$l1)**$scaling)*($hours0*60+$mins0));
  $mins_calc=$mins0;
  } else {
  $nTasks=$l1a[$l1]/$nMols_per_Task;
  $hours_calc=$hours0;
  $mins_calc=$mins0;
  }

  $sTasks=sprintf("%05d",$nTasks);

  $path = "$workDir";
  $path =~ s/\//\\\//g;

  system "sed -e \'s/##$l2l##/$path\\\/templates-liq\\\/$l2s/g\' -e \'s/##$l1l##/$l1s/g\' $workDir/input.template > input.$out_id.inp";
  system "echo Submitting $l2s Molecules at $nCores cores for $hours_calc:$mins_calc";
  
  # Args[default]: queue[-1] cores cores/node[-1] node-type[-1] mem/core[-1] hours mins name inp out #links exe
  system "submit_cp2k.pl -1 $nTasks $nTasks_per_Node -1 -1 $hours_calc $mins_calc SSc-$l2s-$l1s-$nCores input.$out_id.inp output.$out_id.$sTasks.out 1 $cp2k_ver";
  
  chdir("$workDir/$l1d");
 
  } #level2

 chdir("$workDir");

} #level1

