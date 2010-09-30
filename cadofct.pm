#!/usr/bin/perl -w

# Copyright 2008 Pierrick Gaudry, Emmanuel Thome, Paul Zimmermann,
#                Jeremie Detrey
#
# This file is part of CADO-NFS.
#
# CADO-NFS is free software; you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 2.1 of the License, or (at your option)
# any later version.
#
# CADO-NFS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with CADO-NFS; see the file COPYING.  If not, write to the Free
# Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
# 02110-1301, USA.

package cadofct;
use Exporter;
our @ISA= qw(Exporter);
our @EXPORT=qw(%param $tab_level &read_machines &read_param &do_polysel_bench
&do_sieve_bench &do_factbase &do_init &do_task &banner &info &last_line
&format_dhms);

use strict;
use warnings;

use File::Basename;
use Cwd qw(abs_path);
use List::Util qw[min];
use POSIX qw(ceil);
use Math::BigInt;


###############################################################################
# Message and error handling ##################################################
###############################################################################

# Current level of indentation
our $tab_level = 0;

# Should we use colors (for terminal output) or not?
my $use_colors = defined $ENV{CADO_COLOR} ? $ENV{CADO_COLOR} : 1;

# Terminal width
my $term_cols = 80;

# Pads a string with spaces to match the speficied width
sub pad {
    my $str = "" . shift;
    my $w   = shift() - length $str;
    $w = 0 if $w < 0;
    return $str . (" " x $w);
}

# Formats a message by inserting an indented prefix in front of each line
sub format_message {
    my $prefix = ("    " x $tab_level) . shift;
    my $prefix_raw = $prefix;

    $prefix_raw =~ s/\033\[[^m]*m//g;
    $prefix = $prefix_raw unless $use_colors;

    my @msg;
    for (split /\n/, shift) {
        next if /^$/;

#    my $len = $term_cols - length $prefix_raw;
#
#        while (length $_ > $len) {
#            my $i = rindex $_, " ", $len-1;
#            if (($i = rindex $_, " ", $len-1) >= 0 ||
#                ($i =  index $_, " ", $len)   >= 0) {
#                push @msg, (substr $_, 0, $i);
#                $_ = substr $_, $i+1;
#            } else {
#                last;
#            }
#        }
        push @msg, $_;
    }
    s/^/$prefix/, s/$/\n/ for @msg;
    return join "", @msg;
}

my $log_fh;     # filehandle for logging.

# Message function
sub info {
    my $text=shift;
    print STDERR format_message("\033[01;32mInfo\033[01;00m:", $text);
    print $log_fh format_message("Info:", $text) if defined($log_fh);
}

# Banner function
sub banner {
    $tab_level = 0;
    info (("-" x ($term_cols-6))."\n".shift()."\n".("-" x ($term_cols-6)));
}

# Warning hook
$SIG{__WARN__} = sub {
    my $text=shift;
    print $log_fh format_message("Warning:", $text) if defined($log_fh);
    warn         format_message("\033[01;33mWarning\033[01;00m:", $text);
};

# Error hook
$SIG{__DIE__}  = sub {
    die @_ if $^S; 
    my $text=shift;
    print $log_fh format_message("Error:", $text) if defined($log_fh);
    die          format_message("\033[01;31mError\033[01;00m:", $text);
};


###############################################################################
# Parameters ##################################################################
###############################################################################

# Default parameters
# This list gives:
#  - the preferred ordering for parameters;
#  - the default values (if any).
my @default_param = (
    # global
    wdir         => undef,
    bindir      => undef,
    name         => undef,
    machines     => undef,
    n            => undef,
    parallel     => 0,

    # polyselect using Kleinjung (kj) algorithm
    degree       => 5,
    kjkeep       => 100,
    kjkmax       => 10,
    kjincr       => 60,
    kjl          => 7,
    kjM          => 1e25,
    kjpb         => 256,
    kjp0max      => 100000,
    kjadmin      => undef,
    kjadmax      => undef,
    kjadrange    => 1e7,
    kjdelay      => 120,
    selectnice   => 10,

    # sieve
    rlim         => 8000000,
    alim         => 8000000,
    lpbr         => 29,
    lpba         => 29,
    mfbr         => 58,
    mfba         => 58,
    rlambda      => 2.3,
    alambda      => 2.3,
    I            => 13,
    excess       => 1,
    qmin         => 12000000,
    qrange       => 1000000,
    checkrange   => 1000000,

    delay        => 120,
    sievenice    => 19,
    keeprelfiles => 0,
    sieve_max_threads => 2,
    ratq	 => 0,

    # filtering
    keep         => 160, # should be 128+skip
    excesspurge  => 1,
    keeppurge    => 160,
    maxlevel     => 15,
    cwmax        => 200,
    rwmax        => 200,
    ratio        => 1.5,
    bwstrat      => 1,
    skip         => 32,

    # linalg
    linalg       => 'bwc',
    bwmt         => 2,
    mpi		 => 0,
    hosts	 => "",
    bwthreshold  => 64,
    bwtidy       => 1,
    bwc_interval => 1000,
    bwc_mm_impl => 'bucket',
    bwc_interleaving => 0,

    # characters
    nkermax      => 30,
    nchar        => 50,

    # holy grail
    expected_factorization => undef,

    # logfile
    logfile => undef,
);

# Hash for the parameters, global to avoid passing it to each function
our %param = @default_param; # initialize to default values

# Build the ordered list of parameters
my @param_list;

while (@default_param) {
    push @param_list, shift @default_param;
    shift @default_param;
}



# Parses command-line and configuration file parameters.
# The second parameter is a hash of options:
#  - `strict' specifies whether the checking should be strict or not
#             (i.e. die in case of parsing errors)
sub read_param {
    my ($param, $opt) = (shift, shift);
    my $count_args = 0;

    my %files_read = ();

    @_ = map { [$_, 0] } @_;

    ARGS : while (defined ($_ = shift)) {
        die unless ref $_ eq 'ARRAY';
        my $secondary = $_->[1];
        $_=$_->[0];
        next if $files_read{$_};
        $files_read{$_}=1;
        if (-d $_) {
            info "interpreting directory $_ as wdir=$_";
            $_ = "wdir=$_";
        }
        for my $p (@param_list) {
            if (/^$p=(.*)$/) {
                my $v=$1;
                if ($secondary && $p =~ /^(wdir)$/) {
                    if (!defined($param->{$p}) || $param->{$p} ne $v) {
                        warn "$_ ignored in secondary input file\n";
                    }
                    next ARGS;
                }
                $param->{$p} = $v;
                $count_args++;
                my $f;
                if ($p eq 'wdir') {
                    $f = "$1/param";
                } elsif ($p eq 'name' && defined($param->{'wdir'})) {
                    $f = "$param->{'wdir'}/$param->{'name'}.param";
                }
                if (defined($f) && -f $f && !$files_read{$f}) {
                    $count_args--;
                    info "Reading extra parameters from $f\n";
                    unshift @_, [$f,1]
                }
                next ARGS;
            }
        }
        if (/^params?=(.*)$/) {
            warn "Paramfile is not the first argument !\n" if ($count_args);
            my $file = $1;
            open FILE, "< $file"
                or die "Cannot open `$file' for reading: $!.\n";
            my @args;
            while (<FILE>) {
                s/^\s+//; s/\s*(#.*)?$//;
                next if /^$/;
                if (/^(\w+)=(.*)$/) {
                    push @args, "$1=$2";
                    next;
                }
                die "Cannot parse line `$_' in file `$file'.\n"
                    if $opt->{'strict'};
            }
            close FILE;
            unshift @_, map { [$_,$secondary] } @args;
        } elsif (-f $_) {
            unshift @_, [ "param=$_", $secondary ];
        } else {
            die "Unknown argument: `$_'.\n" if $opt->{'strict'};
        }
    }

    # sanity check: old config files may still have true/false values
    while (my ($k, $v) = each %$param) {
        next unless defined $v;
        $param->{$k} = 0 if $v eq 'false';
        $param->{$k} = 1 if $v eq 'true';
    }

    # checking mandatory parameters
    if ($opt->{'strict'}) {
        for my $k ("wdir", "name", "kjadmin", "kjadmax") {
            die "The parameter `$k' is mandatory.\n" if !$param->{$k};
        }
        die "The parameter `machines' is mandatory for parallel mode.\n"
            if $param->{'parallel'} && !$param->{'machines'};

        if (!$param->{'parallel'} && !$param->{'bindir'}) {
            warn "Taking current script's directory as `bindir'.\n";
            $param->{'bindir'} = abs_path(dirname($0));
        }
    }

    # substitute `name' instead of `%s' into `wdir'
    $param->{'wdir'} =~ s/%s/$param->{'name'}/g if ($param->{'wdir'});

    # `prefix' is a shorthand for `$param->{'wdir'}/$param->{'name'}'
    $param->{'prefix'} = "$param->{'wdir'}/$param->{'name'}";
}

# Dumps the list of parameters to a file
sub write_param {
    my ($file) = @_;

    open FILE, "> $file"
        or die "Cannot open `$file' for writing: $!.\n";
    for my $k (@param_list) {
        my $v = $param{$k};
        print FILE "$k=$v\n" if $v;
    }
    close FILE;
}

# Global hash for the machine descriptions
my %machines;

# Reads the machine description file for parallel computing
sub read_machines {
    my $file = shift;
    die "No machine description file was defined.\n" if !$param{'machines'};
    if ( $file ) {
        open FILE, "< $file"
            or die "Cannot open `$param{'machines'}' for reading: $!.\n";
    } else {
        # read from given location if given as an absolute file,
        # otherwise understand as a location relative to the working
        # directory.
        my $m = $param{'machines'};
        if ($m !~ m{/}) {
            info "interpreting filename $m as relative to wdir $param{'wdir'}\n";
            $m = "$param{'wdir'}/$m";
        }
        open FILE, "< $m" or die "Cannot open `$m' for reading: $!.\n";
    }

    my %vars = ();
	
    while (<FILE>) {
        s/^\s+//; s/\s*(#.*)?$//;
        next if /^$/;

        if (/^\[(\S+)\]$/) {
            %vars = ( cluster => $1 );
        } elsif (/^(\w+)=(.*)$/) {
            $vars{$1} = $2;
            $param{'bindir'} = $2 if ($1 eq "bindir");
            if ($1 eq "tmpdir") {
                    my $wdir = $param{'wdir'};
                    $wdir = abs_path(dirname($wdir))."/".basename($wdir);
                    my $tmpdir = abs_path(dirname($2))."/".basename($2);
                    die "tmpdir must be different of wdir in parallel mode.\n"
                            if $wdir eq $tmpdir;
            }
        } elsif (s/^(\S+)\s*//) {
            my $host = $1;
            my %desc = %vars;
            while (s/^(\w+)=(\S*)\s*//) {
                $desc{$1} = $2;
            }
            die "Cannot parse line `$_' in file `$param{'machines'}'.\n"
                if !/^$/;

            for my $k ("tmpdir", "bindir") {
                die "The parameter `$k' is mandatory in $param{'machines'}.\n"
					if !$desc{$k};
            }
            $desc{'tmpdir'} =~ s/%s/$param{'name'}/g;
            $desc{'cores'}  = 1 unless defined $desc{'cores'};
            $desc{'poly_cores'} = $desc{'cores'} 
                    unless defined $desc{'poly_cores'};
            $desc{'mpi'} = 0 unless defined $desc{'mpi'};
            if ( $desc{'mpi'} ) {
                    die "the directory wdir or bindir don't exist on $host.\n"
                    if ( remote_cmd($host, "env test -d $param{'wdir'} && test ".
                                    "-d $param{'bindir'}")->{'status'} );
            }
            while ( $desc{'mpi'} ) {
                    $desc{'mpi'}--;
                    $param{'mpi'}++;
                    $param{'hosts'} .= "$host,";
            }
            $desc{'prefix'} = "$desc{'tmpdir'}/$param{'name'}";
            $desc{'files'}  = {}; # List of files uploaded to the host
            $machines{$host} = \%desc;
        } else {
            die "Cannot parse line `$_' in file `$param{'machines'}'.\n";
        }
    }

    close FILE;

    if ( $param{'mpi'} ) {
            chop $param{'hosts'};
            open FILE, "< $param{'bindir'}/linalg/bwc/bwc.pl"
                    or die "Cannot open `$param{'bindir'}/linalg/bwc/bwc.pl' for reading: $!.\n";
            while (<FILE>) {
                    next unless /^my \$mpiexec='';$/;
                    die "CADO-NFS has not been compiled with MPI flag.\n".
                    "Please add the path of the MPI library in the file local.sh ".
                    "(for example: MPI=/usr/lib64/openmpi) and recompile.\n";
            } 
            close FILE;
    } else {
            if (exists($ENV{'OAR_JOBID'})) {
                    open FILE, "> $param{'machines'}.tmp"
                            or die "Cannot open `$param{'machines'}.tmp' for writing: $!.\n";
                    print FILE "tmpdir=$vars{'tmpdir'}\n";
                    print FILE "bindir=$param{'bindir'}\n";
                    print FILE "mpi=1\n";
                    close FILE;
                    cmd( "uniq $ENV{'OAR_NODEFILE'} >> $param{'machines'}.tmp" );
                    read_machines( "$param{'machines'}.tmp" );
            }
            # TODO: Support other scheduling environments.
    }
}



###############################################################################
# Executing local and remote shell commands ###################################
###############################################################################

# Log file for commands
my $cmdlog;

# Runs the command $cmd.
# The second argument is an optional pointer to a hash table of options:
#  - `log'  specifies that the command should be logged into the command log;
#  - `kill' specifies that if the command fails (non-zero exit status)
#           we should report the error and die immediately.
# The return value is another pointer to a hash table:
#  - `out'    is the captured standard output of the command;
#  - `status' is the exit status of the command.
#  - `stderr' where to redirect stderr
#  - `stdout_and_stderr' redirect both stdout and stderr to given file.
#  - `append_output' whether to append to redirection file in one of the
#     two cases above.
sub cmd {
    my ($cmd, $opt) = @_;

    my $error_file;
    my $redir='>';
    if ($opt->{'append_output'}) {
        $redir='>>'
    }
    if (defined($_=$opt->{'stderr'})) {
        $error_file=$_;
        $cmd .= " 2$redir$_";
    } elsif (defined($_=$opt->{'stdout_stderr'})) {
        $error_file=$_;
        $cmd .= " $redir$_ 2>&1";
    }

    if ($cmdlog && $opt->{'log'}) {
        open LOG, ">> $cmdlog"
            or die "Cannot open `$cmdlog' for writing: $!.\n";
        print LOG "$cmd\n";
        close LOG;
    }

    my $out    = `$cmd`;
    my $status = $? >> 8;

    if ($? && $opt->{'kill'}) {
        my $diagnostic= "Command `$cmd' terminated unexpectedly" .
                " with exit status $status.\n";
        if (!defined($opt->{'stdout_stderr'})) {
            if ($out) {
                $out =~ s/^/STDOUT: /m;
                $diagnostic .= "$out\n";
            } else {
                $diagnostic .= "STDOUT: none.\n";
            }
        }
        if ($error_file) {
            $diagnostic .= `tail $error_file | sed -e 's,^,STDERR: ,'`;
        }
        die $diagnostic;
    }

    return { out => $out, status => $status };
}

# Runs the command $cmd on the remote host $host.
# The second argument is an optional pointer to a hash table of options:
#  - `nolog'   specifies that the command should not be logged into the
#              command log;
#  - `kill'    specifies that if the command fails (non-zero exit status)
#              or times out we should report the error and die immediately;
#  - `timeout' specifies how many seconds to wait before giving up.
# The return value is another pointer to a hash table:
#  - `out'     is the captured standard output of the command;
#  - `status'  is the exit status of the command;
sub remote_cmd {
    my ($host, $cmd, $opt) = @_;

    $opt->{'timeout'} = 30 unless $opt->{'timeout'};

    $cmd =~ s/"/\\"/g;

    # don't ask for a password: we don't want to fall into interactive mode
    # all the time (especially not in the middle of the night!)
    # use public-key authentification instead!
	my $rsh = 'ssh';
    if (exists($ENV{'OAR_JOBID'})) {
        $rsh="/usr/bin/oarsh";
    }

    $cmd = "env $rsh -q ".
           "-o ConnectTimeout=$opt->{'timeout'} ".
           "-o ServerAliveInterval=".int($opt->{'timeout'}/3)." ".
           "-o PasswordAuthentication=no ".
           "$host \"env sh -c '$cmd'\" 2>&1";

    my $ret = cmd($cmd, $opt);

    warn "Remote access to `$host' failed".
         ($ret->{'out'} ? ":\n$ret->{'out'}\n" : ".\n")
        if $ret->{'status'} == 255;

    return $ret;
}

###############################################################################
# Remote jobs #################################################################
###############################################################################

# Reads a job status file
# Format:
# <host> <pid> <threads> <file> <param1> <param2> ...
# The PID is "done" when the job is finished.
sub read_jobs {
    my ($file) = @_;
    my $jobs = [];

    if (!-f $file) {
        info "No job status file found. Creating empty one.\n";
        return $jobs;
    }

    open FILE, "< $file"
        or die "Cannot open `$file' for reading: $!.\n";
    while (<FILE>) {
        s/^\s+//; s/\s*(#.*)?$//;
        next if /^$/;

        if (s/^(\S+)\s+(\d+|done)\s+(\d+)\s+(\S+)\s*//) {
            my @param = split;
            my %job = ( host => $1, threads => $3, file => $4, param => \@param );
            $job{'pid'} = $2 unless $2 eq "done";
            push @$jobs, \%job;
        } else {
            die "Cannot parse line `$_' in file `$file'.\n";
        }
    }
    close FILE;

    return $jobs;
}

# Dumps job status to a file
sub write_jobs {
    my ($jobs, $file) = @_;

    open FILE, "> $file"
        or die "Cannot open `$file' for writing: $!.\n";
    for my $job (@$jobs) {
        print FILE "$job->{'host'} ".($job->{'pid'} ? $job->{'pid'} : "done").
				   " $job->{'threads'} $job->{'file'} ".
				   join(" ", @{$job->{'param'}})."\n";
    }
    close FILE;
}



# Job description string, with padded fields
sub job_string {
    my ($job) = @_;
    my @name = split (/\./, basename($job->{'file'}));
    my $str = "$name[0] ";  
    $str .= pad($job->{'host'}, 16);
    $str .= " ".pad($_, 8) for @{$job->{'param'}};
    return $str;
}

# Checks if a remote job is still alive.
# Returns 1 if the job is alive, 0 if the job is dead, or -1 if the remote
# connection failed (e.g. timeout).
sub is_job_alive {
    my ($job) = @_;

    # Check if there is a running process of that PID.
    my $ret = remote_cmd($job->{'host'}, "env kill -0 $job->{'pid'}");
    return -1 if $ret->{'status'} == 255;
    return 0  if $ret->{'status'};

    # Using lsof, check if this process is accessing the $job->{'file'} file.
    # We need to call readlink here to get the _absolute_ path of the file,
    # as returned by lsof.
    # FIXME: on some hosts (e.g. our trojans) lsof is not in PATH
    if (0) {
      $ret = remote_cmd($job->{'host'}, "env lsof -Fn -a -p$job->{'pid'} -d1 | ".
                        "env grep ^n\\`env readlink -f $job->{'file'}\\`\$");
      return -1 if $ret->{'status'} == 255;
      return 0  if $ret->{'status'};
    }

    return 1;
}

# Gets the status of a job.
# Returns 1 if the job is still running, 0 if the job finished, or -1 if
# the job died.
# The job is said to have finished if the last line of the output file
# matches against a given pattern.
# If we're unable to determine the job status, we asssume that it's
# still running.
sub job_status {
    my ($job, $pattern) = @_;
    my $status;

    my $alive = is_job_alive($job);
    my $ret;
    if ($job->{'file'} =~ /\.gz$/) {
            $ret = remote_cmd($job->{'host'},
                    "env zcat $job->{'file'} | tail -n1 2>&1");
    } else {	
            $ret = remote_cmd($job->{'host'},
                    "env tail -n1 $job->{'file'} 2>&1");
    }

    if ($ret->{'status'} == 255) {
        info "Unknown status. Assuming job is still running.\n";
        $status = 1;
    } elsif (!$ret->{'status'} && $ret->{'out'} =~ /$pattern/) {
        info "Finished!\n";
        $status = 0;
    } elsif (!$ret->{'status'} && $ret->{'out'} =~ /No such file or directory$/) {
        die "The executable was not found. Make sure the `bindir' parameter ".
            "is valid for host `$job->{'host'}'.\n";
    } elsif (!$ret->{'status'} && $ret->{'out'} =~ /BUG/) {
		die $ret->{'out'};
    } else {
        warn "Could not access output file `$job->{'file'}'.\n"
            if $ret->{'status'};

        if ($alive) {
            info $alive == 1 ? "Running...\n" :
                "Unknown status. Assuming job is still running.\n";
            $status = 1;
        } else {
            warn "Dead?!\n";
            $status = -1;
        }
    }

    return $status;
}

# Kills a remote job.
# The $keep argument prevents the output file to be removed on the host.
sub kill_job {
    my ($job, $keep) = @_;
    info "Killing job:  ".job_string($job)."\n";
    $tab_level++;
    if (is_job_alive($job) == 1) {
        remote_cmd($job->{'host'}, "env kill -9 $job->{'pid'}");
        remote_cmd($job->{'host'}, "env rm -f $job->{'file'}") unless $keep;
    }
    $tab_level--;
}



# Sends a file to a remote host.
sub send_file {
    my ($host, $file) = @_;
    my $m = $machines{$host};
    my $ret;

    # If the file is already supposed to be here, just check
    if ($m->{'files'}->{$file}) {
        $ret = remote_cmd($host, "env test -e $m->{'tmpdir'}/$file 2>&1");
        return unless $ret->{'status'};
        delete $m->{'files'}->{$file};
    }

    info "Sending `$file' to `$host'...";
    $tab_level++;

    # Try to upload the file
    $ret = remote_cmd($host, "env mkdir -p $m->{'tmpdir'} 2>&1");
    $ret = cmd("env rsync --timeout=30 $param{'wdir'}/$file ".
               "$host:$m->{'tmpdir'}/ 2>&1", { log => 1 })
        unless $ret->{'status'};

    if ($ret->{'status'}) {
        warn "$ret->{'out'}\n";
    } else {
        $m->{'files'}->{$file} = 1;
    }
    $tab_level--;
}

# Retrieves the output file of a finished job from a remote host.
# The $keep argument prevents the file to be removed on the host.
# Returns 1 if the download was successful, -1 if the file was not there, or
# 0 if another error occurred (meaning that we might want to try again).
sub get_job_output {
    my ($job, $keep) = (@_);

    my $ret = cmd("env rsync --timeout=30 $job->{'host'}:$job->{'file'} ".
                  "$param{'wdir'}/ 2>&1", { log => 1 });
    my $status = 1;
    if ($ret->{'status'}) {
        my @out = split /\n/, $ret->{'out'};
        warn "$out[0]\n";
        $status = $out[0] =~ /No such file or directory/ ? -1 : 0;
    } elsif (!$keep) {
        remote_cmd($job->{'host'}, "env rm -f $job->{'file'}");
    }
    return $status;
}



###############################################################################
# Miscellaneous functions #####################################################
###############################################################################

# Counts the line of a file, _not_ matching a given regexp.
sub count_lines {
    my ($f, $re) = @_;

    my $n = 0;
    if ($f =~ /\.gz$/) {
            $n= cmd ( "zcat $f | grep -v '#' | wc -l" )->{'out'};
            chomp $n;
            return $n;
    }
    # This seems to be a tad faster than grep -v '$re' | wc -l, so...
    open FILE, "< $f"
       	or die "Cannot open `$f' for reading: $!.\n";
    while (<FILE>) {
        $n++ unless $re && /$re/;
    }
    close FILE;

    return $n;
}

# Returns the first line of a file.
sub first_line {
    my ($f) = @_;
    open FILE, "< $f" or die "Cannot open `$f' for reading: $!.\n";
    $_ = <FILE>;
    close FILE;
    chomp;
    return $_;
}

# Returns the last line of a file.
sub last_line {
    my ($f) = @_;
    my $last = "";
    if ($f =~ /\.gz$/) {
            $last= cmd ("zcat $f | tail -n 1" )->{'out'};
            chomp $last;
            return $last;
    }
    open FILE, "< $f" or die "Cannot open `$f' for reading: $!.\n";
	
    # That should be enough to catch the last line
    seek FILE, -512, 2;
    $last = $_ while <FILE>;
    close FILE;
    chomp $last;
    return $last;
}

# This is _ugly_: the siever takes some parameters via the polynomial file.
# The job of this function is to maintain the sieving parameters this
# $name.poly file up to date.
# TODO: Find a cleaner way to do this! (e.g. command-line parameters for las)
sub append_poly_params {
    my @list = qw(rlim alim lpbr lpba mfbr mfba rlambda alambda);
    my $list = join "|", @list;

    # Strip the parameters at the end of the poly file, in case they
    # have changed
    open IN, "< $param{'prefix'}.poly"
        or die "Cannot open `$param{'prefix'}.poly' for reading: $!.\n";
    open OUT, "> $param{'prefix'}.poly_tmp"
        or die "Cannot open `$param{'prefix'}.poly_tmp' for writing: $!.\n";
    while (<IN>) {
        print OUT "$_" unless /^($list):\s*/;
    }
    close IN;

    # Append the parameters to the poly file
    print OUT "$_: $param{$_}\n" for @list;
    close OUT;

    cmd("env mv -f $param{'prefix'}.poly_tmp $param{'prefix'}.poly",
        { kill => 1 });
}

sub local_time {
    my $job= shift;
    $cmdlog = "$param{'prefix'}.cmd";
    open LOG, ">> $cmdlog" or die "Cannot open `$cmdlog' for writing: $!.\n";
    print LOG "# Starting $job on " . localtime() . "\n";
    close LOG;
}
    
sub format_dhms {
    my $sec = shift;
    my ($d, $h, $m);
    $d = int ( $sec / 86400 ); $sec = $sec % 86400;
    $h = int ($sec / 3600 ); $sec = $sec % 3600;
    $m = int ($sec / 60 ); $sec = $sec % 60;
    return "$d"."d:$h"."h:$m"."m:$sec"."s"; 
}
    
###############################################################################
# Distributed tasks ###########################################################
###############################################################################

# Scans a list of ranges and merges overlapping ones.
sub merge_ranges {
    my ($ranges) = @_;
    my @merged = ();

    for my $r (sort { $a->[0] <=> $b->[0] } @$ranges) {
        my ($a, $b) = @$r;
        if (!@merged || $a > $merged[-1]->[1]) {
            push @merged, $r;
            next;
        } elsif ($b > $merged[-1]->[1]) {
            $merged[-1]->[1] = $b;
        }
    }

    return \@merged;
}

# Finds a hole of a bounded size in a given interval, excluding already
# listed ranges. Returns () if no suitable hole was found.
sub find_hole {
    my ($min, $max, $len, $ranges) = @_;

    die "Invalid range: `$min-$max'.\n" if $max && $min >= $max;

    # Remove ranges lying completely before [$min,$max]
    shift @$ranges while scalar @$ranges && $ranges->[0]->[1] < $min;

    # Insert dummy [$min,$min] range if there is room at the beginning
    unshift @$ranges, [$min,$min]
        unless scalar @$ranges && $ranges->[0]->[0] <= $min;

    # The hole starts right after the first range
    # We allocate a full $len-sized hole first
    my $a = $ranges->[0]->[1];
    my $b = $a + $len;

    # Truncate the hole if needed
    $b = $max              if $max                && $max              < $b;
    $b = $ranges->[1]->[0] if scalar @$ranges > 1 && $ranges->[1]->[0] < $b;

    # Make sure the hole is a proper range
    return ($a, $b) if $a < $b;
    return ();
}



# This function is the common factor of the polynomial selection and sieving
# codes.
# Its only argument is a huge hash with keys:
#  - `task'       is the name of the task, like "polysel" or "sieve", used to
#                 name job status files;
#  - `title'      is the title of the task, to be displayed in a flashy banner;
#  - `suffix'     is the common suffix of the job output files of the task
#                 (i.e. "kjout" for polynomial selection, or "rels" for
#                 sieving);
#  - `extra'      is an optional suffix, to match extra files when recovering
#                 job output files (i.e. "freerels" for sieving);
#  - `files'      is a list of the files that need to be sent to each host;
#  - `pattern'    is a regexp to match the last line of a completed job output
#                 file;
#  - `min', `max' the bounds on the range to process; `max' is facultative
#                 (meaning that the range will grow until we have enough data);
#  - `len'        the maximal size of a range to be processed by a job;
#  - `partial'    is a flag specifying if we can import partial job output
#                 files: if a job died, can we still use its output?
#  - `keep'       is a flag specifying if we should leave the job output files
#                 on the hosts;
#  - `delay'      is the number of seconds to wait between each polling of the
#                 job status;
#  - `check'      is a subroutine which checks the integrity of a job output
#                 file; it should remove the file if it's invalid, and fix it
#                 if it's not complete; this function takes a second parameter
#                 specifiying if the check should be exhaustive or not;
#  - `progress'   is a subroutine which will print the progress of the current
#                 task;
#  - `is_done'    is a subroutine which checks if we're done or not; it takes
#                 the list of ranges;
#  - `cmd'        is a subroutine which, given a range and a host machine
#                 description, returns the command to run the task on the host.
sub distribute_task {
    my ($opt) = @_;

    banner $opt->{'title'};
    local_time $opt->{'title'};
	$opt->{'gzip'}=0 if (! $opt->{'gzip'});

    # Make sure that all the output files that are already here are correct
    opendir DIR, $param{'wdir'}
        or die "Cannot open directory `$param{'wdir'}': $!\n";
    my $suffix = $opt->{'suffix'}.'\.[\de.]+-[\de.]+(|\.gz)';
    $suffix .= "|$opt->{'extra'}" if $opt->{'extra'};
    my @files = grep /^$param{'name'}\.($suffix)$/,
                     readdir DIR;
    closedir DIR;

    if (@files) {
        info "Checking previous files...\n";
        $tab_level++;
        # We don't do exhaustive checking here, it's too slow...
        # We assume the files are here for a good reason!
        &{$opt->{'check'}}($_, 0) for (map "$param{'wdir'}/$_", sort @files);
        $tab_level--;
    }



    while (1) {
        my $jobs = [];

        # See what's already running, and retrieve possibly finished data
        if ($param{'parallel'}) {
            $jobs = read_jobs("$param{'prefix'}.$opt->{'task'}_jobs");

            my @running = grep  defined $_->{'pid'}, @$jobs;
            my @done    = grep !defined $_->{'pid'}, @$jobs;
            my @new_jobs;

            # Check the status of all running jobs
            if (@running) {
                info "Checking all running jobs...\n";
                $tab_level++;
                for my $job (@running) {
                    info "Checking job: ".job_string($job)."\n";
                    $tab_level++;
                    my $status = job_status($job, $opt->{'pattern'});
                    if ($status == 1) {
                        # Job is still alive: keep it in the list
                        push @new_jobs, $job;
                    } elsif ($status == 0 || $opt->{'partial'}) {
                        # Job is (partially) terminated: mark it as done
                        delete $job->{'pid'};
                        push @done, $job;
                    } else {
                        # Job is dead: remove its output file on the host
                        remote_cmd($job->{'host'}, "env rm -f $job->{'file'}");
                    }
                    $tab_level--;
                }
                $tab_level--;
            }

            # Retrieve files of finished jobs
            if (@done) {
                info "Retrieving job data...\n";
                $tab_level++;
                for my $job (@done) {
                    info "Retrieving `".basename($job->{'file'})."' ".
                         "from `$job->{'host'}'...\n";
                    $tab_level++;

                    my $file = "$param{'wdir'}/".basename($job->{'file'});
                    if (-f $file) {
                        warn "`$file' already exists. ".
                             "Assuming it is the same.\n";
                    } else {
                        my $status = get_job_output($job, $opt->{'keep'});
                        if ($status == 1) {
                            # Output file was downloaded: exhaustive check
                            &{$opt->{'check'}}($file, 1);
                        } elsif ($status == 0) {
                            # Can't get output file: let's try again next time
                            push @new_jobs, $job;
                        } else {
                            # File is not there: too bad...
                        }
                    }
                    $tab_level--;
                }
                $tab_level--;
            }

            $jobs = \@new_jobs;
            write_jobs($jobs, "$param{'prefix'}.$opt->{'task'}_jobs");
        }

        # Scan ranges
        my $ranges = [];

        # First, scan the job output files
        opendir DIR, $param{'wdir'}
            or die "Cannot open directory `$param{'wdir'}': $!\n";
        my @files = grep /^$param{'name'}\.$opt->{'suffix'}\.[\de.]+-[\de.]+(|\.gz)$/,
                         readdir DIR;
        closedir DIR;
        push @$ranges, map { /\.([\de.]+)-([\de.]+)(|\.gz)$/; [$1, $2] } @files;
        $ranges = merge_ranges($ranges);

        # Keep a copy for later
        my $file_ranges = [ map [@$_], @$ranges ];

        # Add the ranges from running or done jobs
        my $good_jobs = [ grep { $_->{'file'} =~ /$param{'name'}\./ } @$jobs ];
        push @$ranges, map { my @p = @{$_->{'param'}}; \@p } @$good_jobs;
        $ranges = merge_ranges($ranges);

        # Start new job(s) (parallel mode)
        if ($param{'parallel'}) {
            info "Starting new jobs...\n";
            $tab_level++;

            HOST : for my $h (keys %machines) {
                my $m = $machines{$h};
                my $cores= $m->{'cores'};
                $cores = $m->{'poly_cores'} if ($opt->{'task'} eq "polysel"); 
				
                # How many free cores on this host?
                my $busy_cores = 0;
                foreach (@$jobs) {
                    $busy_cores += $_->{'threads'} if $_->{'host'} eq $h;
                }
                my $n = $cores - $busy_cores;
                next if $n < 1;

                # Send files and skip to next host if not all files are here
                # (don't do this as an anonymous loop, as I've seen it
                # behave oddly).
                for my $f (@{$opt->{'files'}}) {
                    send_file($h, $f);
                }
                for my $f (@{$opt->{'files'}}) {
                    next if $m->{'files'}->{$f};
                    warn "$h does not have file $f, skipping host\n";
                    next HOST;
                }

                my $nth = $opt->{'max_threads'};
                while ($n > 0) {
                    $n -= $opt->{'max_threads'};
                    $nth = $n + $opt->{'max_threads'} if $n < 0;
                    my @r = find_hole($opt->{'min'}, $opt->{'max'},
                                      $opt->{'len'}, $ranges);

                    # No hole was found. But maybe we are waiting for another
                    # job to finish, on a host which is unreachable...
                    # So instead of staying idle, let's be redundant!
                    # (This patch is sponsored by your local energy provider!)
                    #if (!@r) {
                    #    $ranges = [ map [@$_], @$file_ranges ];
                    #    @r = find_hole($opt->{'min'}, $opt->{'max'},
                    #                   $opt->{'len'}, $ranges);
                    #}

                    # Still no hole? Well then, we're truly finished!
                    last HOST unless @r;

                    my $job = { host  => $h,
                            threads => $nth,
                            file  => "$m->{'prefix'}.$opt->{'suffix'}.".
                                         "$r[0]-$r[1]",
                            param => \@r };
                    $job->{'file'} .= ".gz" if $opt->{'gzip'};

                    info "Starting job: ".job_string($job)."\n";
                    $tab_level++;
                    my $cmd = &{$opt->{'cmd'}}(@r, $m, $nth, $opt->{'gzip'}).
                            " & echo \\\$!";
                    my $ret = remote_cmd($h, $cmd, { log => 1 });
                    if (!$ret->{'status'}) {
                        chomp $ret->{'out'};
                        $job->{'pid'} = $ret->{'out'};
                        push @$jobs, $job;
                        push @$ranges, \@r;
                        $ranges = merge_ranges($ranges);
                    }
                    $tab_level--;
                }
            }

            write_jobs($jobs, "$param{'prefix'}.$opt->{'task'}_jobs");
            $tab_level--;
        }

        # Print the progress of the task
        &{$opt->{'progress'}}($file_ranges) if $opt->{'progress'};

        # This might be enough to exit the loop now
        last if &{$opt->{is_done}}($file_ranges);



        # Start new job (sequential mode)
        if (!$param{'parallel'} && (my @r = find_hole($opt->{'min'}, $opt->{'max'},
                                $opt->{'len'}, $ranges)))
        {
                # XXX What is bwmt doing in here ???
                my $mt = $param{'bwmt'};
                $mt=$1*$2 if ($mt =~ /^(\d+)x(\d+)$/);
                my $nth = min ( $opt->{'max_threads'}, $mt );

                info "Starting job: ".pad($r[0], 8)." ".pad($r[1], 8)."\n";
                $tab_level++;
                my $cmd = &{$opt->{'cmd'}}(@r, $machines{'localhost'}, $nth, $opt->{'gzip'});
                cmd($cmd, { log => 1, kill => 1 });
                my $check_cmd = "$param{'prefix'}.$opt->{'suffix'}.$r[0]-$r[1]";
                $check_cmd .= ".gz" if $opt->{'gzip'};
                &{$opt->{'check'}}($check_cmd, 1); # Exhaustive testing!
                $tab_level--;
        }

        # Wait for a bit before having another go
        if ($param{'parallel'}) {
            info "Waiting for $opt->{'delay'} seconds before ".
                 "checking again...\n";
            sleep $opt->{'delay'};
        }
    }

    # A bit of cleaning on slaves
    if (! $opt->{'bench'}) {
        info "Cleaning up...\n";
        $tab_level++;
        # Kill jobs
        my $jobs = read_jobs("$param{'prefix'}.$opt->{'task'}_jobs");
        for my $job (@$jobs) {
            kill_job($job, $opt->{'partial'});
            next unless $opt->{'partial'};
            $tab_level++;
            my $file = "$param{'wdir'}/".basename($job->{'file'});
            if (-f $file) {
                warn "`$file' already exists. Assuming it is the same.\n";
            } else {
                get_job_output($job, $opt->{'keep'});
                &{$opt->{'check'}}("$param{'wdir'}/".basename($job->{'file'}), 1);
                # TODO: For now, the extra relations are imported back in the
                # working directory, but they are not used. Feeding them to
                # duplicates/singleton would take some more time and we don't
                # really care about them since we already have enough
                # relations.
            }
            $tab_level--;
        }
        unlink "$param{'prefix'}.$opt->{'task'}_jobs";
        # Remove files
        while (my ($h, $m) = each %machines) {
            my $files = join " ", (map "$m->{'tmpdir'}/$_", @{$opt->{'files'}});
            remote_cmd($h, "env rm -f $files");
        }
        $tab_level--;
    }
}



###############################################################################
# Tasks #######################################################################
###############################################################################

# List of tasks with their dependencies:
#  - `name'   is the task name;
#  - `dep'    is the list of tasks on which the current task depends:
#             if one of these tasks is more recent, we also have to
#             reschedule the current task;
#  - `req'    is the list of order-only dependencies: we have to complete
#             all these tasks before scheduling the current task (but
#             there is no notion of "more recent" here);
#  - `param'  is the list of parameters on which the current task depends:
#             reschedule the task is a parameter has changed;
#  - `files'  is the list of suffix patterns of the files generated by this
#             task, used for automatic cleanup;
#  - `resume' specifies that the task can be resumed if all its dependencies
#             are up to date;
#  - `dist'   specifies that the task can be distributed, used to kill all
#             the jobs when cleaning up.
#
# Some fields will be added during the execution of the script:
#  - `rdep'    is the list of reverse dependecies, i.e. the tasks that depend
#              on this one;
#  - `rreq'    is the list of reverse order-only dependecies;
#  - `visited' is used by graph traversal algorithms;
#  - `done'    is the time at which the task has been completed (if any).
my %tasks = (
    init      => { },

    polysel   => { name   => "polynomial selection",
                   dep    => ['init'],
                   param  => ['degree', 'kjM', 'kjl', 'kjkeep', 'kjkmax',
                              'kjincr', 'kjpb', 'kjp0max', 'kjadmin',
                              'kjadmax'],
                   files  => ['kjout\.[\de.]+-[\de.]+', 'poly', 'poly_tmp'],
                   resume => 1,
                   dist   => 1 },

    factbase  => { name   => "factor base",
                   dep    => ['polysel'],
                   param  => ['alim'],
                   files  => ['roots', 'makefb\.stderr'] },

    freerels  => { dep    => ['factbase'],
                   files  => ['freerels.gz', 'freerel\.stderr'] },

    sieve     => { name   => "sieve and purge",
                   dep    => ['polysel'],
                   req    => ['factbase', 'freerels'],
                   param  => ['excess'],
                   files  => ['rels\.[\de.]+-[\de.]+(|\.gz)', 'rels\.tmp',
                              'nodup\.gz', 'dup1\.stderr', 'dup2\.stderr',
                              'purged', 'purge\.stderr'],
                   resume => 1,
                   dist   => 1 },

    merge     => { name   => "merge",
                   dep    => ['sieve'],
                   param  => ['keep', 'maxlevel', 'cwmax', 'rwmax',
                              'ratio', 'bwstrat'],
                   files  => ['merge\.his', 'merge\.stderr'] },

    # replay shouldn't appear as a step in its own right. It's a bug.
    replay    => { name   => "replay",
                   dep    => ['merge'],
                   files  => ['index', 'small.bin', 'replay\.stderr'],
                   param  => ['skip'], },

    linalg    => { name   => "linear algebra",
                   dep    => ['replay'],
                   param  => [ qw/bwmt bwthreshold linalg
                               bwc_interval
                               bwc_mm_impl
                               bwc_interleaving/],
                   files  => ['bwc', 'bwc\.stderr', 'bl', 'bl\.stderr',
							  'W\d+'] },

    chars     => { name   => "characters",
                   dep    => ['linalg'],
                   param  => ['nchar'],
                   files  => ['ker', 'characters\.stderr'] },

    sqrt	  => { name	  => "square root",
				   dep    => ['chars'],
                   param  => ['nkermax'],
                   files  => ['dep\.\d+', 'dep\.alg\.\d+', 'dep\.rat\.\d+',
                              'sqrt\.stderr', 'fact\.\d+',
							  'fact', 'allfactors'] }
);

# Initialize empty arrays
for my $v (values %tasks) {
    for (qw(dep req rdep rreq param files)) {
        $v->{$_} = [] unless defined $v->{$_};
    }
}

# Build reverse dependencies
while (my ($k, $v) = each %tasks) {
    push @{$tasks{$_}->{'rdep'}}, $k for @{$v->{'dep'}};
    push @{$tasks{$_}->{'rreq'}}, $k for @{$v->{'req'}};
}



# Runs a task, after possibly running the tasks on which it depends first.
sub do_task {
    my ($t) = @_;
    my $task = $tasks{$t};

    # Do nothing if the task was already completed
    return if $task->{'done'};

    # First, do all tasks on which this one depends
    do_task($_) for (@{$task->{'dep'}}, @{$task->{'req'}});

    # Call the corresponding do_* function
    # (we need to allow symbolic refs for that)
    {
        no strict 'refs';
        &{"do_$t"}();
    }

    # Put a timestamp file
    open FILE, "> $param{'prefix'}.${t}_done"
        or die "Cannot open `$param{'prefix'}.${t}_done' for writing: $!.\n";
    close FILE;
    $task->{'done'} = (stat("$param{'prefix'}.${t}_done"))[9] # modificaton time
        or die "Cannot stat `$param{'prefix'}.${t}_done': $!.\n";
}



###############################################################################
# Initialization ##############################################################
###############################################################################

sub do_init {
    banner "Initialization";

    # Getting configuration
    info "Reading the parameters...\n";
    $tab_level++;
    read_param(\%param, { strict => 1 }, @ARGV);
    $tab_level--;

    if ($param{'parallel'}) {
        info "Reading the machine description file...\n";
        $tab_level++;
        read_machines();
        $tab_level--;
    } else {
        $machines{'localhost'} = { tmpdir  => $param{'wdir'},
                                 bindir => $param{'bindir'},
                                 prefix  => $param{'prefix'} };
    }

    info "Initializing the working directory...\n";
    $tab_level++;
    # Create working directory if not there
    cmd("env mkdir -p $param{'wdir'} 2>&1", { kill => 1 })
        if !-d $param{'wdir'};
    $tab_level--;

    if (defined($param{'logfile'})) {
        open $log_fh, ">$param{'logfile'}" or die "$param{'logfile'}: $!";
    }

    # Check if there is already some stuff relative to $name in $wdir
    # First thing is $name.n. If it is not there, we consider that
    # everything is obsolete, anyway.
    my $recover = 0;
    if (-f "$param{'prefix'}.n") {
        info "There is already some data relative to `$param{'name'}' ".
             "in the working directory. Trying to recover...\n";
        $tab_level++;
        $recover = 1;

        open FILE, "< $param{'prefix'}.n"
            or die "Cannot open `$param{'prefix'}.n' for reading: $!.\n";
        $_ = <FILE>;
        close FILE;
        chomp;
        die "Cannot parse `$param{'prefix'}.n'.\n" unless /^n:\s*(\d+)$/;

        if (!$param{'n'}) {
            $param{'n'} = $1;
        } elsif ($param{'n'} != $1) {
            warn "The contents of `$param{'name'}.n' are inconsistent ".
                 "with the given parameter `n'. Aborting recovery.\n";
            $recover = 0;
        }

        $tab_level--;
    }

    # If something was done here before, retrieve the parameters to see
    # from where we should start again
    my %param_diff;
    if ($recover && -f "$param{'prefix'}.param") {
        eval {
            my %param_old;
            read_param(\%param_old, { strict => 0 },
                       "param=$param{'prefix'}.param");
            for (keys %param) {
            		$param_diff{$_} =$param{$_} ne $param_old{$_}
						if (exists($param_old{$_}));
			}
        };
    }

    if (!$recover) {
        # Read n if not given on command line
        if (!$param{'n'}) {
            info "The parameter `n' was not specified. Please enter the ".
                 "number to factor:\n";
            $param{'n'} = <STDIN>;
            chomp $param{'n'};
        }

        # Create $name.n in $wdir
        open FILE, "> $param{'prefix'}.n"
            or die "Cannot open `$param{'prefix'}.n' for writing: $!.\n";
        print FILE "n: $param{'n'}\n";
        close FILE;
    }

    local_time(basename($0));

    # Timestamp the task with the date of the last modification of $name.n
    $tasks{'init'}->{'done'} = (stat("$param{'prefix'}.n"))[9] # modification time
        or die "Cannot stat `$param{'prefix'}.n': $!.\n";



    # Task recovery using the dependency graph
    # Topological sort and traversal of the task dependency graph, to check
    # which tasks are up to date.
    my @queue = ("init");
    my @cleanup;
    while (my $t = shift @queue) {
        my $task = $tasks{$t};

        # Skip if already visited
        next if $task->{'visited'};

        # First, check that all previous nodes have been visited, otherwise
        # skip this node for now
        for (map $tasks{$_}->{'visited'}, (@{$task->{'dep'}}, @{$task->{'req'}})) {
            next unless $_;
        }

        # Visit this node and push the next ones in the queue
        $task->{'visited'} = 1;
        push @queue, (@{$task->{'rdep'}}, @{$task->{'rreq'}});

        my $done   = $task->{'done'};
        my $resume = $task->{'resume'} && $recover;

        # Is there already a $name.${t}_done file? If so, it must mean that
        # the task has already been done
        if (!$done && -f "$param{'prefix'}.${t}_done") {
            $done = (stat("$param{'prefix'}.${t}_done"))[9] # modification time
                or die "Cannot stat `$param{'prefix'}.${t}_done': $!.\n";
        }

        # Check dependencies
        if ($done || $resume) {
            for (map $tasks{$_}, @{$task->{'dep'}}) {
                if (!$_->{'done'}) {
                    info "$_->{'name'} not flagged as done, flagging ${t} as ".
                      "not done\n"
						if $_->{'name'};
                    undef $done;
                    undef $resume;
                    last;
                }
                if ($done && $_->{'done'} > $done) {
                    info "$_->{'name'}_done newer than ${t}_done\n"
						if $_->{'name'};
                    undef $done;
                    undef $resume;
                    last;
                }
            }
        }

        # Check parameter changes
        if ($done || $resume) {
            for (@{$task->{'param'}}) {
                if ($param_diff{$_}) {
                    info "Parameters changed for ${t}\n";
                    undef $done;
                    undef $resume;
                    last;
                }
            }
        }

        # If the task is up to date or can be resumed, we're done for now
        if ($done) {
            info "Nothing to be done for $task->{'name'}.\n" if $task->{'name'};
            $task->{'done'} = $done;
        } else {
            delete $task->{'done'};
        }
        next if $done || $resume;

        # Otherwise, add to the clean-up list
        my $files = join "|", (@{$task->{'files'}}, "${t}_done");
        opendir DIR, $param{'wdir'}
            or die "Cannot open directory `$param{'wdir'}': $!\n";
        my @files = grep /^$param{'name'}\.($files)$/,
                         readdir DIR;
        closedir DIR;
        push @cleanup, { task => $t, files => \@files } if @files;
    }

    # Clear the `visited' field of each node
    delete $_->{'visited'} for values %tasks;

    # Cleaning up everything in one go.
    if (@cleanup) {
        # Make sure that the user is ready to cleanup!
        my $list = join "; ", (grep { defined $_ }
                               (map $tasks{$_->{'task'}}->{'name'}, @cleanup));
        $list =~ s/^(.*);/$1; and/;
        warn "I will clean up the following tasks: $list.\n";
        warn "Are you OK to continue? [y/l/N] (30s timeout)\n";
        warn "(l: recover linear algebra with checkpoint, clean the other tasks)\n";
        my $r = "";
        eval {
            local $SIG{'ALRM'} = sub { die "alarm\n" }; # NB: \n required
            alarm 30;
            $r = <STDIN>;
            alarm 0;
        };
        if ($@) {
            die unless $@ eq "alarm\n"; # propagate unexpected errors
        }
        chomp $r;
        die "Aborting...\n" unless $r =~ /^(y|l)/i;

        for (@cleanup) {
            my $t     = $_->{'task'};
            my $files = $_->{'files'};
            my $task  = $tasks{$t};

            # Clean up old files...
            if ( $task->{'name'} ) {
                if ( $task->{'name'} eq "linear algebra" ) {
                    next if $r eq "l";
                }
            }
            info "Cleaning up $task->{'name'}..." if $task->{'name'};
            $tab_level++;
            for (map "$param{'wdir'}/$_", sort @$files) {
                unlink $_ if -f;
                cmd("env rm -rf $_") if -d;
            }
            $tab_level--;

            # ... and kill old jobs
            if ($task->{'dist'} && -f "$param{'prefix'}.${t}_jobs") {
                info "Killing old $task->{'name'} jobs..."
                    if $task->{'name'} && -s "$param{'prefix'}.${t}_jobs";
                $tab_level++;
                my $jobs = read_jobs("$param{'prefix'}.${t}_jobs");
                kill_job($_) for @$jobs;
                unlink "$param{'prefix'}.${t}_jobs";
                $tab_level--;
            }
        }
    }

    # Dump parameters into $name.param
    write_param("$param{'prefix'}.param");

    # Update parameters in the $name.poly file if needed
    append_poly_params() if $tasks{'polysel'}->{'done'};
}



###############################################################################
# Polynomial selection ########################################################
###############################################################################

my $polysel_check = sub {
    my ($f) = @_;
    if (! -f $f) {
        warn "File `$f' not found.\n";
        return;
    }

    my %poly;
    open FILE, "< $f"
        or die "Cannot open `$f' for reading: $!.\n";
    while (<FILE>) {
        if (/^No polynomial found/) {
            warn "No polynomial in file `$f'.\n".
			     "check [kj]M value.\n";
            close FILE;
            return;
        }
        $poly{$1} = $2 if /^(\w+):\s*([\w\-.]+)$/;
    }
    close FILE;
        
    # Remove invalid files
    for (qw(n skew Y1 Y0), map "c$_", (0 .. $param{'degree'})) {
        if (!defined $poly{$_}) {
            warn "File `$f' is incomplete (missing `$_'). Removing...\n";
            unlink $f;
            return;
        }
    }
    if ($poly{'n'} != $param{'n'}) {
        warn "File `$f' is invalid (different `n'). Removing...\n";
        unlink $f;
    }
};

my $polysel_cmd = sub {
    my ($a, $b, $m, $max_threads, $gzip) = @_;
    return "env nice -$param{'selectnice'} ".
           "$m->{'bindir'}/polyselect/polyselect ".
           "-keep $param{'kjkeep'} ".
           "-kmax $param{'kjkmax'} ".
           "-incr $param{'kjincr'} ".
           "-l $param{'kjl'} ".
           "-M $param{'kjM'} ".
           "-pb $param{'kjpb'} ".
           "-p0max $param{kjp0max} ".
           "-admin $a ".
           "-admax $b ".
           "-degree $param{'degree'} ".
           "< $m->{'prefix'}.n ".
           "> $m->{'prefix'}.kjout.$a-$b ".
           "2>&1";
};

sub do_polysel {
    my $polysel_is_done = sub {
        my ($ranges) = @_;
        for (@$ranges) {
            next     if $_->[1] <  $param{'kjadmax'};
            last     if $_->[0] >  $param{'kjadmin'};
            return 1 if $_->[0] <= $param{'kjadmin'} &&
                        $_->[1] >= $param{'kjadmax'};
        }
        return 0;
    };

    my $polysel_progress = sub {
        my ($ranges) = @_;
        my ($min, $max) = ($param{'kjadmin'}, $param{'kjadmax'});

        my $total = 0;
        for (@$ranges) {
            my @r = ($_->[0] < $min ? $min : $_->[0],
                $_->[1] > $max ? $max : $_->[1]);
            $total += $r[1] - $r[0] if $r[0] < $r[1];
        }
        $total = (100 * $total) / ($max - $min);

        info "Total interval coverage: ".sprintf("%3.0f", $total)." %.\n";
    };

    distribute_task({ task     => "polysel",
                      title    => "Polynomial selection",
                      suffix   => "kjout",
                      files    => ["$param{'name'}.n"],
                      pattern  => '^(# generated|No polynomial found)',
                      min      => $param{'kjadmin'},
                      max      => $param{'kjadmax'},
                      len      => $param{'kjadrange'},
                      delay    => $param{'kjdelay'},
                      check    => $polysel_check,
                      progress => $polysel_progress,
                      is_done  => $polysel_is_done,
                      cmd      => $polysel_cmd,
					  max_threads => 1 });

    info "All done!\n";

    # Choose best according to the Murphy value
    my $Emax;
    my $best;

    opendir DIR, $param{'wdir'}
        or die "Cannot open directory `$param{'wdir'}': $!\n";
    my @files = grep /\.kjout\.[\de.]+-[\de.]+$/,
                     readdir DIR;
    closedir DIR;

    for my $f (map "$param{'wdir'}/$_", sort @files) {
        open FILE, "< $f"
            or die "Cannot open `$f' for reading: $!.\n";
        my $last;
        my $line;
        while ($line=<FILE>) {
            if ($line =~ /Murphy/ ) {
                $last = $line;
                last;
            }
        }
        close FILE;

        next unless $last && $last =~ /\)=(.+)$/;
        if (!defined $Emax || $1 > $Emax) {
            $Emax = $1;
            $best = $f;
        }
    }

    die "No polynomial was found in the given range!\n".
        "Please increase the range or the [kj]M value.\n"
      unless defined $Emax;

    # Copy the best polynomial
    info "The best polynomial is from `".basename($best)."' (E = $Emax).\n";
    $tab_level++;
    cmd("env cp -f $best $param{'prefix'}.poly 2>&1",
        { log => 1, kill => 1 });
    $tab_level--;

    # Append sieving parameters to the poly file
    open FILE, ">> $param{'prefix'}.poly"
        or die "Cannot open `$param{'prefix'}.poly' for writing: $!.\n";
    print FILE "$_: $param{$_}\n"
        for qw(rlim alim lpbr lpba mfbr mfba rlambda alambda);
    close FILE;
}

sub do_polysel_bench {
	my $last = shift;

    my $polysel_is_done = sub {
        my ($ranges) = @_;
        my ($min, $max) = ($param{'kjadmin'}, $param{'kjadmax'});

        my $total = 0;
        for (@$ranges) {
            my @r = ($_->[0] < $min ? $min : $_->[0],
                     $_->[1] > $max ? $max : $_->[1]);
            $total += $r[1] - $r[0] if $r[0] < $r[1];
        }
        $total = ceil ($total / $param{'kjadrange'});
		my $total_cores=0;
		foreach (keys %machines) {
			$total_cores += $machines{$_}{'poly_cores'};
		}
        my $size = count_lines("$param{'prefix'}.polysel_jobs", "$param{'name'}\.");
		my $total_jobs = ceil (($max-$min)/$param{'kjadrange'});
		if ($last) {
			return 1 if $total >= $total_jobs + $size;
		} else {
			return 1 if $total > $total_jobs - $total_cores + $size;
		}
		return 0;
    };

    distribute_task({ task     => "polysel",
                      title    => "Polynomial selection",
                      suffix   => "kjout",
                      files    => ["$param{'name'}.n"],
                      pattern  => '^(# generated|No polynomial found)',
                      min      => $param{'kjadmin'},
                      max      => $param{'kjadmax'},
                      len      => $param{'kjadrange'},
                      delay    => $param{'kjdelay'},
                      check    => $polysel_check,
                      is_done  => $polysel_is_done,
                      cmd      => $polysel_cmd,
                      bench    => 1,
                      max_threads => 1 });

    if ($last) {
    	info "All done!\n";
    } else {
    	info "Switch to next configuration...\n";
    }
}



###############################################################################
# Factor base #################################################################
###############################################################################

sub do_factbase {
    info "Generating factor base...\n";
    $tab_level++;

    my $cmd = "$param{'bindir'}/sieve/makefb ".
              "-poly $param{'prefix'}.poly ".
              "> $param{'prefix'}.roots ";
    cmd($cmd, { log => 1, kill => 1,
            stderr=>"$param{'prefix'}.makefb.stderr" });
    $tab_level--;
}



###############################################################################
# Free relations ##############################################################
###############################################################################

sub do_freerels {
    info "Computing free relations...\n";
    $tab_level++;

    my $cmd = "$param{'bindir'}/sieve/freerel ".
              "-poly $param{'prefix'}.poly ".
              "-fb $param{'prefix'}.roots ".
              "> $param{'prefix'}.freerels ";

    cmd($cmd, { log => 1, kill => 1,
            stderr=>"$param{'prefix'}.freerel.stderr" });
	cmd("gzip $param{'prefix'}.freerels");
    $tab_level--;
}



###############################################################################
# Sieve and purge #############################################################
###############################################################################

my $sieve_cmd = sub {
    my ($a, $b, $m, $max_threads, $gzip) = @_;
    my $cmd = "env nice -$param{'sievenice'} ".
        "$m->{'bindir'}/sieve/las ".
        "-I $param{'I'} ".
        "-poly $m->{'prefix'}.poly ".
        "-fb $m->{'prefix'}.roots ".
        "-q0 $a ".
        "-q1 $b ".
        "-mt $max_threads ";
    $cmd .= "-ratq " if ($param{'ratq'});
    $cmd .=	"-out $m->{'prefix'}.rels.$a-$b";
    $cmd .= ".gz" if ($gzip);
    $cmd .= " > /dev/null 2>&1";
    return $cmd;
};

sub do_sieve {
    my $nrels      = 0;
    my $last_check = 0;

    my $import_rels = sub {
        my ($f) = @_;
        my $n = count_lines($f, '^#');
        $nrels += $n;
        info "Imported $n relations from `".basename($f)."'.\n";
    };

    # XXX. No. choose a way -- separate packages, whatever. But an
    # anonymous function of this size is a no-go.
    my $sieve_check = sub {
        my ($f, $full) = @_;

        unless (-f $f) {
            warn "File `$f' not found, check not done.\n";
            return;
        }


        return &$import_rels($f) if $f =~ /\.freerels.gz$/;
        my $is_gzip;
        $is_gzip=1 if $f =~ /\.gz$/;

        my $check = $f;
        if (!$full) {
            $check = "$param{'prefix'}.rels.tmp";
            # Put the first 10 relations into a temp file
            if ($is_gzip) {
                open FILE, "zcat $f|"
                    or die "Cannot open `$f' for reading: $!.\n";
            } else {
                open FILE, "< $f"
                    or die "Cannot open `$f' for reading: $!.\n";
            }
            open TMP, "> $check"
                or die "Cannot open `$check' for writing: $!.\n";
            my $n = 10;
            while (<FILE>) {
                $n--, print TMP $_ unless /^#/;
                last unless $n;
            }
            close FILE;
            close TMP;
        }

        # Check relations
        my $ret = cmd("$param{'bindir'}/utils/check_rels ".
                      "-poly $param{'prefix'}.poly $check > /dev/null 2>&1");
        unlink "$check" unless $full;

        # Remove invalid files
        if ($ret->{'status'} == 1) {
            my $msg="File `$f' is invalid (check_rels failed).";
            if ($ENV{'CADO_DEBUG'}) {
                my $nf = "$f.error";
                $msg .= " Moving to $nf\n";
                warn $msg;
                rename $f, $nf;
            } else {
                $msg .= " Deleting.\n";
                warn $msg;
                unlink $f;
            }
            close FILE;
            return;
        } elsif ($ret->{'status'}) {
            # Non-zero, but not 1? Something's wrong, bail out
            die "check_rels exited with unknown error code ", 
                 $ret->{'status'}, ", aborting."
        }

        # If this is a partial (i.e. incomplete) file, we need to adjust
        # the range of covered special q's
        if (last_line($f) !~ /^# (Total \d+ reports|Warning: truncated)/) {
            if ($is_gzip) {
                cmd ("gzip -d $f", { kill => 1});
                basename($f) =~ /^$param{'name'}\.rels\.([\de.]+)-([\de.]+)\.gz$/;
                $f = "$param{'prefix'}.rels.$1-$2";
            }
            open FILE, "+< $f"
                or die "Cannot open `$f' for update: $!.\n";

            # TODO: Since the file is truncated, we assume that the last
            # reported special q was not completely sieved, so we remove it.
            # Maybe we can still save it, but is it worth the trouble?
            my @lastq;
            my $pos = 0;
            while (<FILE>) {
                # Keep track of the last two special q's
                if (/^### q=(\d+): roots?/) {
                    shift @lastq if scalar @lastq == 2;
                    push @lastq, { q => $1, pos => $pos };
                }
                $pos = tell FILE;
            }

            # Less than two special q's in this file: nothing to recover
            if (scalar @lastq < 2) {
                warn "File `$f' contains no useable data. Deleting...\n";
                close FILE;
                unlink $f;
                return;
            }

            # Truncate the file and add a marker at the end
            truncate FILE, $lastq[-1]->{'pos'};
            seek FILE, $lastq[-1]->{'pos'}, 0;
            print FILE "# Warning: truncated file\n";
            close FILE;

            # Rename the file to account for the truncated range
            basename($f) =~ /^$param{'name'}\.rels\.([\de.]+)-([\de.]+)$/;
            my @r = ($1, $lastq[0]->{'q'}+1);
            info "Truncating `".basename($f)."' to range $r[0]-$r[1]...\n";
            $tab_level++;
            cmd("env mv -f $f $param{'prefix'}.rels.$r[0]-$r[1]", { kill => 1 });
            $f = "$param{'prefix'}.rels.$r[0]-$r[1]";
            if ($is_gzip) {
                cmd ("gzip $f", { kill => 1});
                $f .= ".gz";
            }
            $tab_level--;
        }

        # The file is clean: we can import the relations now
        &$import_rels($f);
    };



    my $sieve_progress = sub {
        info "Running total: $nrels relations.\n";
    };



    my $sieve_is_done = sub {
        # Check only every $param{'checkrange'} relations
        return 0 if $nrels - $last_check < $param{'checkrange'};
        $last_check = $nrels;

        # Get the list of relation files
        opendir DIR, $param{'wdir'}
            or die "Cannot open directory `$param{'wdir'}': $!\n";

        my $pat=qr/^$param{'name'}\.(rels\.[\de.]+-[\de.]+|freerels)\.gz$/;

        my @files = grep /$pat/, readdir DIR;
        closedir DIR;
        mkdir "$param{'prefix'}.nodup"
        unless (-d "$param{'prefix'}.nodup");
        my $nslices = 4;
        for (my $i=0; $i < $nslices; $i++) {
                mkdir "$param{'prefix'}.nodup/$i"
                unless (-d "$param{'prefix'}.nodup/$i");
        }
        opendir DIR, "$param{'prefix'}.nodup/0/"
            or die "Cannot open directory `$param{'prefix'}.nodup/0/': $!\n";
        my @old_files = grep /$pat/, readdir DIR;
        closedir DIR;
        my %old_files;		
        $old_files{$_} = 1 for (@old_files);
        my @new_files;
        for (@files) {
                push @new_files, $_ unless (exists ($old_files{$_}));
        }

        # print number of primes in factor base
        if (scalar @files >= 2) {
            my $f = $files[0];
            $f = $files[1]
                if $files[0] =~ /^$param{'name'}\.freerels.gz$/;
            $f = "$param{'wdir'}/".$f;
            open FILE, "zcat $f|"
                or die "Cannot open `$f' for reading: $!.\n";
            my $i=0;
            while (<FILE>) {
                if ( $_ =~ /^# (Number of primes in \S+ factor base = \d+)$/ ) {
                    info "$1\n";
                    $i++;
                    last if $i==2; 
                }
            }
            close FILE;
        }

        banner "Duplicate and singleton removal";
        # Remove duplicates
        info "Removing duplicates...";
        $tab_level++;
        if (@new_files) {
                my $new_files = join " ",
                (map "$param{'wdir'}/$_", sort @new_files);
                info "split new files in $nslices slices...";
                cmd("$param{'bindir'}/filter/dup1 ".
                        "-out $param{'prefix'}.nodup $new_files ",
                        { log => 1, kill => 1,
                                stdout_stderr=>"$param{'prefix'}.dup1.stderr" });
        }
        {
            my $name="$param{'prefix'}.subdirlist";
            open FILE, "> $name" or die "$name: $!";
            print FILE join("\n", map { "$param{'name'}.nodup/$_"; } (0..$nslices-1));
            close FILE;
        }
        {
            # Put basenames of relation files in list.
            my $name="$param{'prefix'}.filelist";
            open FILE, "> $name" or die "$name: $!";
            print FILE join("\n", map { m{([^/]+)$}; $1; } @files);
            close FILE;
        }

        my $n = 0;
        my @allfiles;
        my $K = int ( 100 + (1.2 * $nrels / $nslices) );
        for (my $i=0; $i < $nslices; $i++) {
                info "removing duplicates on slice $i...";
                cmd("$param{'bindir'}/filter/dup2 ".
                        "-K $K -out $param{'prefix'}.nodup/$i ".
                        "-filelist $param{'prefix'}.filelist ".
                        "-basepath $param{'prefix'}.nodup/$i ",
                        { log => 1, kill => 1,
                                stdout_stderr => "$param{'prefix'}.dup2.stderr",
                        });
                my $f = "$param{'prefix'}.dup2.stderr";
                open FILE, "< $f"
                        or die "Cannot open `$f' for reading: $!.\n";
                while (<FILE>) {
                        if ( $_ =~ /^\s+(\d+) remaining relations/ ) {
                                # shiftwidth 8 in perl scripts is
                                # totally nuts.
                                $n += $1;
                                last;
                        }
                }
                close FILE;
        }

        info "Number of relations left: $n.\n";
        $tab_level--;

        # Remove singletons
        info "Removing singletons...";
        $tab_level++;
        my $ret = cmd("$param{'bindir'}/filter/purge ".
                "-poly $param{'prefix'}.poly -keep $param{'keeppurge'} ".
                "-excess $param{'excesspurge'} ".
                "-nrels $n -out $param{'prefix'}.purged ".
                "-basepath $param{'wdir'} " .
                "-subdirlist $param{'prefix'}.subdirlist ".
                "-filelist $param{'prefix'}.filelist ",
                { log => 1,
                        stdout_stderr => "$param{'prefix'}.purge.stderr"
                });
        if ($ret->{'status'}) {
                info "Not enough relations! Continuing sieving...\n";
                $tab_level--;
                return 0;
        }

        # Get the number of rows and columns from the .purged file
        my ($nrows, $ncols) = split / /, first_line("$param{'prefix'}.purged");
        my $excess = $nrows - $ncols;
        if ($excess < $param{'excess'}) {
            info "Not enough relations! Continuing sieving...\n";
            $tab_level--;
            return 0;
        }
		
        info "Nrows: $nrows; Ncols: $ncols; Excess: $excess.\n";
        $tab_level--;
        # note: the nodup.gz file is no longer used now.
        # Thus we may spare the hassle of creating it.
        info "Join all no duplicate files into one file...";
        cmd("cat $param{'prefix'}.purgefiles | xargs zcat ".
            "| gzip --best > $param{'prefix'}.nodup.gz ",
            { log => 1, kill => 1 });
        $tab_level++;
        info "clean directory nodup...";
        cmd("rm -rf $param{'prefix'}.nodup");
        $tab_level--;

        return 1;
    };
    
    distribute_task({ task     => "sieve",
                      title    => "Sieve",
                      suffix   => "rels",
                      extra    => "freerels.gz",
                      gzip     => 1,
                      files    => ["$param{'name'}.poly",
                                   "$param{'name'}.roots"],
                      pattern  => '^# Total \d+ reports',
                      min      => $param{'qmin'},
                      len      => $param{'qrange'},
                      partial  => 1,
                      keep     => $param{'keeprelfiles'},
                      delay    => $param{'delay'},
                      check    => $sieve_check,
                      progress => $sieve_progress,
                      is_done  => $sieve_is_done,
                      cmd      => $sieve_cmd,
                      max_threads => $param{'sieve_max_threads'} });

    info "All done!\n";
}


sub do_sieve_bench {
    my $max_rels = shift;
    my $last = shift;
    my $nrels      = 0;
    my $max_files;

    my $import_rels = sub {
        my ($f) = @_;
        my $n = count_lines($f, '^#');
        $nrels += $$max_rels[1] * $n / $param{'qrange'};
        info "Imported $n relations from `".basename($f)."'.\n" if $n > 0;
    };

    my $sieve_check = sub {
        my ($f, $full) = @_;

        unless (-f $f) {
                warn "File `$f' not found, check not done.\n";
                return;
        }

        my $is_gzip;
        $is_gzip=1 if $f =~ /\.gz$/;
        my $check = $f;
        if (!$full) {
            $check = "$param{'prefix'}.rels.tmp";
            # Put the first 10 relations into a temp file
            if ($is_gzip) {
                open FILE, "zcat $f|"
                    or die "Cannot open `$f' for reading: $!.\n";
            } else {
                open FILE, "< $f"
                    or die "Cannot open `$f' for reading: $!.\n";
            }
            open TMP, "> $check"
                or die "Cannot open `$check' for writing: $!.\n";
            my $n = 10;
            while (<FILE>) {
                $n--, print TMP $_ unless /^#/;
                last unless $n;
            }
            close FILE;
            close TMP;
        }

        # Check relations
        my $ret = cmd("$param{'bindir'}/utils/check_rels ".
                      "-poly $param{'prefix'}.poly $check > /dev/null 2>&1");
        unlink "$check" unless $full;

        # Remove invalid files
        if ($ret->{'status'} == 1) {
            my $msg="File `$f' is invalid (check_rels failed).";
            if ($ENV{'CADO_DEBUG'}) {
                my $nf = "$f.error";
                $msg .= " Moving to $nf\n";
                warn $msg;
                rename $f, $nf;
            } else {
                $msg .= " Deleting.\n";
                warn $msg;
                unlink $f;
            }
            close FILE;
            return;
        } elsif ($ret->{'status'}) {
            # Non-zero, but not 1? Something's wrong, bail out
            die "check_rels exited with unknown error code ", 
                 $ret->{'status'}, ", aborting."
        }
        # The file is clean: we can import the relations now
        &$import_rels($f);
    };

    my $sieve_progress = sub {
        info "Estimate relations: $nrels.\n";
    };
	
    my $sieve_is_done = sub {
        return 0 if $nrels < $$max_rels[2];

        opendir DIR, $param{'wdir'}
            or die "Cannot open directory `$param{'wdir'}': $!\n";
        my @files = grep /^$param{'name'}\.rels\.[\de.]+-[\de.]+\.gz$/,
                         readdir DIR;
        closedir DIR;
        @files = map { /\.([\de.]+)-[\de.]+\.gz$/; $1 }
        @files;
        @files = sort ( {$a <=> $b} @files );
        $max_files = $files[-1] unless ($max_files);
        my $number_files_total = ( $max_files - $param{'qmin'} ) /
        $$max_rels[1] + 1;
        my $number_files = 1;
        while ($files[0] != $max_files) {
            $number_files++;
            shift @files;
        }
        if ( $number_files == $number_files_total ) {
            return 1;
        } else {
            return 0;
        }
    };

    distribute_task({ task     => "sieve",
                      title    => "Sieve",
                      suffix   => "rels",
                      extra    => "freerels.gz",
                      gzip	   => 1,
                      files    => ["$param{'name'}.poly",
                                   "$param{'name'}.roots"],
                      pattern  => '^# Total \d+ reports',
                      min      => $param{'qmin'},
                      len      => $param{'qrange'},
                      keep     => $param{'keeprelfiles'},
                      delay    => $param{'delay'},
                      check    => $sieve_check,
                      progress => $sieve_progress,
                      is_done  => $sieve_is_done,
                      cmd      => $sieve_cmd,
					  max_threads => $param{'sieve_max_threads'} });

    if ($last) {
    	info "All done!\n";
    } else {
    	info "Switch to next configuration...\n";
    }
}



###############################################################################
# Merge #######################################################################
###############################################################################

my $bwcostmin;

sub do_merge {
    banner "Merge";
    info "Merging relations...\n";
    $tab_level++;

    my $cmd = "$param{'bindir'}/filter/merge ".
              "-out $param{'prefix'}.merge.his ".
              "-mat $param{'prefix'}.purged ".
              "-forbw $param{'bwstrat'} ".
              "-keep $param{'keep'} ".
              "-maxlevel $param{'maxlevel'} ".
              "-cwmax $param{'cwmax'} ".
              "-rwmax $param{'rwmax'} ".
              "-ratio $param{'ratio'} ";

    cmd($cmd, { log => 1, kill => 1, stdout_stderr =>
            "$param{'prefix'}.merge.stderr" });

    if (last_line("$param{'prefix'}.merge.his") =~ /^BWCOSTMIN: (\d+)/) {
        $bwcostmin = $1;
        info "Minimal bwcost: $bwcostmin.\n";
    }

    $tab_level--;
}



###############################################################################
# Replay ######################################################################
###############################################################################

sub do_replay {
    info "Replaying merge history...\n";
    $tab_level++;

    if (!defined $bwcostmin &&
        last_line("$param{'prefix'}.merge.his") =~ /^BWCOSTMIN: (\d+)/) {
        $bwcostmin = $1;
    }

    my $cmd = "$param{'bindir'}/filter/replay ".
              "--binary " .
              "-skip $param{'skip'} " .
              "-his $param{'prefix'}.merge.his ".
              "-index $param{'prefix'}.index ".
              "-purged $param{'prefix'}.purged ".
              "-out $param{'prefix'}.small.bin ".
              (defined $bwcostmin ? "-costmin $bwcostmin " : "");

    cmd($cmd, { log => 1, kill => 1,
              stdout_stderr=>"$param{'prefix'}.replay.stderr "
        });

    my ($nrows, $ncols, $weight);
    open FILE, "< $param{'prefix'}.replay.stderr"
        or die "Cannot open `$param{'prefix'}.replay.stderr' for reading: $!.\n";
    while (<FILE>) {
        $nrows = $1, $ncols = $2 if /^small_nrows=(\d+) small_ncols=(\d+)/;
        $weight = $1             if /^# Weight\(M_small\) = (\d+)/;
    }
    close FILE;
    info "Nrows: $nrows; Ncols: $ncols; Weight: $weight.\n";

    $tab_level--;
}



###############################################################################
# Linear algebra ##############################################################
###############################################################################

sub do_linalg {
    banner "Linear algebra";
    local_time "Linear algebra";

    my $cmd;
    if ($param{'linalg'} eq "bw") {
        die "Old code no longer supported";
    } elsif ($param{'linalg'} eq "bwc") {
        info "Calling Block-Wiedemann (new code)...\n";
        $tab_level++;
        my $mt = $param{'bwmt'};
        if ($mt =~ /^(\d+)$/) {
            $mt = "${mt}x1";
        }

        my $bwc_script = "$param{'bindir'}/linalg/bwc/bwc.pl";

        # Note: $param{'bindir'} is not expanded yet. So if we get it as
        # a variable from the mach_desc file, it won't do. It's better to
        # pass it as a command-line argument to bwc.pl
        
        my $bwc_bindir = "$param{'bindir'}/linalg/bwc";

        # XXX NOTE: This is a despair-mode fallback. It's really not
        # guaranteed to work, even though it's the way I'm sometimes
        # using the script. The ``official'' way is to use the script
        # which is in the build dir, because that one has the @xxx@ stuff
        # replaced (and provides in particular the advantage that
        # bwc_bindir does not need to be specified).
        if (!-x $bwc_script) {
            $bwc_script=abs_path(dirname($0)) . "/linalg/bwc/bwc.pl";
        }
        if (!-x $bwc_script) {
            die "script bwc.pl not found";
        }

        $cmd = "$bwc_script ".
               ":complete " .
               "seed=1 ". # For debugging purposes, we use a deterministic BW
               "thr=$mt ";
        if ( $param{'mpi'} > 1 ) {
            my $a = int ( sqrt($param{'mpi'}) );
            $a-- while ( $param{'mpi'} % $a != 0);
            my $b = $param{'mpi'} / $a;				
            $cmd .= "mpi=$b"."x$a hosts=$param{'hosts'} ";
            # TODO: Support other scheduling environments.
            # TODO: Support non-openmpi command lines.
            if (exists($ENV{'OAR_JOBID'})) {
                $cmd .= 
                "mpi_extra_args='--mca btl_tcp_if_exclude lo,virbr0 --mca plm_rsh_agent oarsh' ";
            } else {
                $cmd .= "mpi_extra_args='--mca btl_tcp_if_exclude lo,virbr0' ";
            }
        } else {
            $cmd .= "mpi=1x1 ";
        }
        $cmd .= "matrix=$param{'prefix'}.small.bin " .
               "nullspace=left " .
               "mm_impl=$param{'bwc_mm_impl'} ".
               "interleaving=$param{'bwc_interleaving'} ".
               "interval=$param{'bwc_interval'} ".
               "mode=u64 mn=64 splits=0,64 ys=0..64 ".
               "wdir=$param{'prefix'}.bwc " .
               "bwc_bindir=$bwc_bindir ";
        cmd($cmd, { log => 1, kill => 1,
                append_output=>1,
                stdout_stderr=>"$param{'prefix'}.bwc.stderr" });

    } elsif ($param{'linalg'} eq "bl") {
        die "No longer supported";
    } else {
        die "Value `$param{'linalg'}' is unknown for parameter `linalg'\n";
    }

    $tab_level--;
}



###############################################################################
# Characters ##################################################################
###############################################################################

my $ndep;

sub do_chars {
    info "Adding characters...\n";
    $tab_level++;

    my $cmd = "$param{'bindir'}/linalg/characters ".
              "-poly $param{'prefix'}.poly ".
              "-purged $param{'prefix'}.purged ".
              "-index $param{'prefix'}.index ".
              "-heavyblock $param{'prefix'}.small.dense.bin ".
              "-nchar $param{'nchar'} ".
              "-out $param{'prefix'}.ker ";

    opendir BWC, "$param{'prefix'}.bwc";
    my @kers = grep { /^K.\d+$/ } readdir BWC;
    closedir BWC;
    if (!scalar @kers) {
        die "No kernel files out of bwc ???";
    }
    $cmd .= " $param{'prefix'}.bwc/$_" foreach @kers;

    cmd($cmd, { log => 1, kill => 1,
            stderr=>"$param{'prefix'}.characters.stderr" });

    $ndep = `awk '/^Wrote/ { print \$2; }' $param{'prefix'}.characters.stderr`;
    chomp($ndep);
    info "$ndep vectors remaining after characters.\n";

    $tab_level--;
}



###############################################################################
# Square root #################################################################
###############################################################################

sub is_prime {
    my $n = shift;
    my $z=0+cmd("$param{'bindir'}/utils/gmp_prob_prime $n")->{'out'};
    return $z;
}

sub primetest_print {
    my $n=shift;
    if (is_prime($n)) {
        return "$n [prime]";
    } else {
        return "$n [composite]";
    }
}

sub do_sqrt {
    banner "Square root";
    local_time "Square root";
    if (!defined($ndep)) {
        $ndep = `awk '/^Wrote/ { print \$2; }' $param{'prefix'}.characters.stderr`;
        chomp($ndep);
    }
    if (!defined($ndep) || $ndep > $param{'nkermax'}) {
        $ndep = $param{'nkermax'};
    }

    # We don't use bigints as hash keys.
    my @prime_factors=();
    my %composite_factors=($param{'n'}=>1);

    {
        # First prepare all deps files
        info "Preparing $ndep dependency files\n";
        my $cmd = "$param{'bindir'}/sqrt/sqrt ".
            "-poly $param{'prefix'}.poly ".
            "-prefix $param{'prefix'}.dep " .
            "-ab " .
            "-purged $param{'prefix'}.purged ".
            "-index $param{'prefix'}.index ".
            "-ker $param{'prefix'}.ker ";

        cmd($cmd, { log => 1, kill => 1,
                append_output=>1,
                stdout_stderr=>"$param{'prefix'}.sqrt.stderr"});
    }

    # later processing does not need re-generation of the .dep files.
    for (my $numdep=0; $numdep<$ndep; $numdep++) {
        my $znumdep=sprintf('%03d', $numdep);
        my $f="$param{'prefix'}.fact.$znumdep";
        info "Testing dependency number $numdep...\n";
        $tab_level++;
        my $cmd = "$param{'bindir'}/sqrt/sqrt ".
            "-poly $param{'prefix'}.poly ".
            "-prefix $param{'prefix'}.dep " .
            "-dep $numdep " .
            "-rat -alg -gcd " .
            "-purged $param{'prefix'}.purged ".
            "-index $param{'prefix'}.index ".
            "-ker $param{'prefix'}.ker ".
            "> $f";

        cmd($cmd, { log => 1, kill => 1,
                append_output=>1,
                stderr=>"$param{'prefix'}.sqrt.stderr"});

        do { $tab_level--; next; } if first_line($f) =~ /^Failed/;

        info "Factorization was successful!\n";
        # only informational.
        cmd("env cp -f $f $param{'prefix'}.fact", { kill => 1 });

        my @factors_thisdep=();
        open FILE, "< $f" 
            or die "Cannot open `$f' for reading: $!.\n";
        while (<FILE>) {
            chomp($_);
            push @factors_thisdep, $_;
        }
        close FILE;

        for my $p (@factors_thisdep) {
            info(primetest_print($p) ."\n");
        }
            
        info "Doing gcds with previously known composite factors\n";

        my @kcomp = keys %composite_factors;
        for my $m (@kcomp) {
            my @nontriv=();
            my $zm = Math::BigInt->new($m);
            for my $p (@factors_thisdep) {
                my $zp = Math::BigInt->new($p);
                my $za = Math::BigInt::bgcd($zp, $zm);
                next if $za->is_one();
                next if $za->bcmp($zm) == 0;
                # We have a non-trivial factor, thus a split of m.
                push @nontriv, $za->bstr();
            }
            if (@nontriv) {
                delete $composite_factors{$m};
                for my $a (@nontriv) {
                    info "non-trivial factor: ".primetest_print($a)."\n";
                    if (is_prime($a)) {
                        push @prime_factors, $a;
                    } else {
                        $composite_factors{$a}=1;
                    }
                }
            }
        }

        my $np = scalar @prime_factors;
        my $nc = scalar keys %composite_factors;
        info "Now: $np prime factors, $nc composite factors\n";

        if ($nc == 0) {
            info "Factorization complete\n";
            $tab_level--;
            last;
        }
        $tab_level--;
    }

    die "No square root was found.\n" unless @prime_factors;

    my $f1 = "$param{'prefix'}.allfactors";
    open FILE, "> $f1" or die "Cannot open `$f1' for writing: $!.\n";
    print FILE "$_\n" for @prime_factors;
    close FILE;
}

close $log_fh if $log_fh;
1;

# vim: set tabstop=8 shiftwidth=4 sta et:
