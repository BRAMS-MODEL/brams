/*
                   University of Illinois/NCSA
                       Open Source License

        Copyright(C) 1998-2006, The Board of Trustees of the
            University of Illinois.  All rights reserved.

                          Developed by:             

                      The PerfSuite Project
         National Center for Supercomputing Applications
            University of Illinois at Urbana-Champaign

                 http://perfsuite.ncsa.uiuc.edu/

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal with the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

+ Redistributions of source code must retain the above copyright notice, 
  this list of conditions and the following disclaimers.
+ Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimers in
  the documentation and/or other materials provided with the distribution.
+ Neither the names of The PerfSuite Project, NCSA/University of Illinois
  at Urbana-Champaign, nor the names of its contributors may be used to
  endorse or promote products derived from this Software without specific
  prior written permission.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
DEALINGS WITH THE SOFTWARE.
*/


#define MAXCOMMSIZ 16

typedef struct {
    pid_t pid;               /* Process ID */

    char comm[MAXCOMMSIZ];   /* The filename of the executable, in
                              * parentheses.  This is visible     
                              * whether or not the executable is  
                              * swapped out.
                              */

    char state;              /* One character from the string
                              * "RSDZTW" where R is running, S is
                              * sleeping in an interruptible wait, D
                              * is sleeping in an uninterruptible
                              * disk sleep, Z is zombie, T
                              * is traced or stopped (on a signal),
                              * and W is paging.
                              */

    pid_t ppid;              /* The PID of the parent. */

    pid_t pgrp;              /* The process group ID of the process. */

    pid_t session;           /* The session ID of the process. */

    int tty;                 /* The tty the process uses. */

    int tpgid;               /* The process group ID of the process
                              * which currently owns the tty that
                              * the process is connected to.
                              */

    unsigned long flags;     /* The flags of the process.
                              * The math bit is decimal 4, and the
                              * traced bit is a decimal 10.
                              */

    unsigned long minflt;    /* The number of minor faults the
                              * process has made, those which have
                              * not required loading a memory page
                              * from disk.
                              */

    unsigned long cminflt;   /* The number of minor faults that the
                              * process and its children have made.
                              */

    unsigned long majflt;    /* The number of major faults the
                              * process has made, those which have
                              * required loading a memory page from
                              * disk.
                              */

    unsigned long cmajflt;   /* The number of major faults that the
                              * process and its children have made.
                              */

    clock_t utime;           /* The number of jiffies that this
                              * process has been scheduled in user
                              * mode.
                              */

    clock_t stime;           /* The number of jiffies that this
                              * process has been scheduled in kernel
                              * mode.
                              */

    clock_t cutime;          /* The number of jiffies that this
                              * process and its children have been
                              * scheduled in kernel mode.
                              */

    clock_t cstime;          /* The number of jiffies that this
                              * process and its children have been
                              * scheduled in kernel mode.
                              */


    long priority;           /* The standard nice value, plus
                              * fifteen.  The value is never
                              * negative in the kernel.
                              */

    long nice;               /* The nice value ranges from 19 (nicest)
                              * to -19 (not nice to others).
                              */

    /* @@ constant 0UL in array.c here... */

    unsigned long itrealvalue;/* The time (in jiffies) before the
                              * next SIGALRM is sent to the process
                              * due to an interval timer.
                              */

    unsigned long starttime; /* Time the process started in jiffies
                              * after system boot.
                              */

    unsigned long vsize;     /* Virtual memory size in bytes */

    unsigned long rss;       /* Resident Set Size: number of pages
                              * the process has in real memory,
                              * minus 3 for administrative purposes.
                              * This is just the pages which count
                              * towards text, data, or stack space.
                              * This does not include pages which
                              * have not been demand-loaded in, or
                              * which are swapped out.
                              */

    unsigned long rlim;      /* Current limit in bytes on the rss of
                              * the process (usually 4,294,967,295).
                              */

    unsigned long startcode; /* The address above which program text
                              * can run.
                              */

    unsigned long endcode;   /* The address below which program text
                              * can run.
                              */

    unsigned long startstack;/* The address of the start of the
                              * stack.
                              */

    unsigned long kstkesp;   /* The current value of esp (stack
                              * pointer), as found in the
                              * kernel stack page for the process.
                              */

    unsigned long kstkeip;   /* The current EIP (instruction
                              * pointer).
                              */

    int signal;              /* The bitmap of pending signals
                              * (usually 0).
                              */

    int blocked;             /* The bitmap of blocked signals
                              * (usually 0, 2 for shells).
                              */

    int sigignore;           /* The bitmap of ignored signals. */

    int sigcatch;            /* The bitmap of catched signals. */

    unsigned long wchan;     /* This is the "channel" in which the
                              * process is waiting.  This is the
                              * address of a system call, and can be
                              * looked up in a namelist if you need
                              * a textual name. (If you have an up-
                              * to-date /etc/psdatabase, then try ps -l
                              * to see the WCHAN field in action)
                              */

    unsigned long nswap;     /* Number of pages swapped - not maintained. */

    unsigned long cnswap;    /* Cumulative nswap for child processes. */

    int exitsignal;          /* Signal to be sent to parent when we die. */

    int processor;           /* Processor number last executed on. */
} ps_procstat_t;

