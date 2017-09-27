#!/usr/bin/env perl
package FEBA_utils;
use strict;
use warnings;

my $_IS_JOB_COMPLETE = "isjobcomplete";

sub wait {
    my $jobId = shift;
    my $sleepTime = shift || 150;

    my $wait = 1;
    while ($wait) {
        sleep($sleepTime);
        my $exitCode = system("$_IS_JOB_COMPLETE $jobId > /dev/null");
        $wait = 0 if !$exitCode;
    }
}

sub create_job_script {
    my $fh = shift;
    my $pe_slots = shift;
    my $ram = shift;
    my $time = shift;
    my $project = shift;
    my $cmds = shift;
    my $modules = shift;
    print $fh "#!/bin/bash -l\n".
              "#\n".
              "#\$ -w e\n".
              "#\$ -cwd\n".
              "#\$ -pe pe_slots $pe_slots\n".
              "#\$ -l ram.c=$ram\n".
              "#\$ -l h_rt=$time\n".
              "#\$ -P $project.p\n".
              "\n";
    if ($modules){
        foreach my $mod (@{$modules}){
            print $fh "module load $mod\n";
        }
        print $fh "\n";
    }
    foreach my $cmd (@{$cmds}) {
        print $fh "$cmd\n";
    }
}

sub submit_job {
    my $qsub_script = shift;
    my $job_name = shift;
    my $debug = shift;
    my $hold_jid = shift;
    
    my $qsub_cmd = "qsub -N $job_name ";
    if ($hold_jid){
        $qsub_cmd .= "-hold_jid $hold_jid ";
    }
    $qsub_cmd .= "$qsub_script ";
    if ($debug){ 
        print STDERR "$qsub_cmd\n";
        return 0;
    }

    my $qsub_out = `$qsub_cmd`;
    if ($qsub_out =~ /Your job (\d+)/){
        return $1;
    } else {
        print "qsub of $qsub_script failed. See below:\n$qsub_out";
    }
}

sub create_job_array_script {
    my $fh = shift;
    my $pe_slots = shift;
    my $ram = shift;
    my $time = shift;
    my $project = shift;
    my $cmds = shift;
    my $modules = shift;
    print $fh "#!/bin/bash -l\n".
              "#\n".
              "#\$ -w e\n".
              "#\$ -cwd\n".
              "#\$ -pe pe_slots $pe_slots\n".
              "#\$ -l ram.c=$ram\n".
              "#\$ -l h_rt=$time\n".
              "#\$ -P $project.p\n".
              "\n".
              "PART_NUM=`printf \"%.2d\" \$((\$SGE_TASK_ID-1))`\n".
              "\n";
    if ($modules){
        foreach my $mod (@{$modules}){
            print $fh "module load $mod\n";
        }
        print $fh "\n";
    }
    foreach my $cmd (@{$cmds}) {
        $cmd =~ s/SUB/\$PART_NUM/g;
        print $fh "$cmd\n";
    }
}

sub submit_job_array {
    my $qsub_script = shift;
    my $job_name = shift;
    my $array_start = shift;
    my $array_end = shift;
    my $array_step = shift;
    my $debug = shift;
    my $hold_jid = shift;
    
    my $qsub_cmd = "qsub -N $job_name -t $array_start-$array_end:$array_step ";
    if ($hold_jid){
        $qsub_cmd .= "-hold_jid $hold_jid ";
    }
    $qsub_cmd .= "$qsub_script ";
    if ($debug){ 
        print STDERR "$qsub_cmd\n";
        return 0;
    }

    my $qsub_out = `$qsub_cmd`;
    if ($qsub_out =~ /Your job-array (\d+).*/){
        return $1;
    } else {
        print "qsub of $qsub_script failed. See below:\n$qsub_out";
    }
}

