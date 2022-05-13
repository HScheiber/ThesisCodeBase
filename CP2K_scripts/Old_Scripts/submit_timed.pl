#!/usr/bin/perl

use POSIX qw/floor/;

######################
$cp2k_ver=0;
$hours0=00;
$mins0=05;
$scaling=1.0;

@l1i=(0,1,2,3,4,5,6);
@l2i=(1);
@l3i=(0);

@l1a=(64,128,256,512,1024,2048,4096);
@l2a=('F','T');
@l3a=('TZV2P-GTH','QZV3P-GTH');

$l1l='NMOL';
$l2l='ALMO';
$l3l='BASS';

$nMols_per_Task=-1; # -1 to fix the number of cores
$nCores=256;
$nTasks_per_Node=16;
$out_id="00";
# LOG: 00: 512 cores
#      01: 1024 cores
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

   $path = "$workDir";
   $path =~ s/\//\\\//g;

   system "sed -e \'s/##$l1l##/$path\\\/templates-liq\\\/$l1s/g\' -e \'s/##$l2l##/$l2s/g\' -e \'s/##$l3l##/$l3s/g\' $workDir/input.template > input.$out_id.inp";

   # Args[default]: queue[-1] cores cores/node[-1] node-type[-1] mem/core[-1] hours mins name inp out #links exe
   system "submit_cp2k_new.pl -1 $nTasks $nTasks_per_Node sandybridge -1 $hours_calc $mins_calc b${l1}${l2}${l3}-$sTasks input.$out_id.inp output.$out_id.$sTasks.out 1 $cp2k_ver";
  
   chdir("$workDir/$l1d/$l2d");
 
  } #level3
 
  chdir("$workDir/$l1d");
 
 } #level2

 chdir("$workDir");

} #level1

