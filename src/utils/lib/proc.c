#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include "procfs.h"

#define MEGABYTE 1048576

#define get_procfs get_procfs_

static int pgsize;

int ps_procstat(ps_procstat_t *);

#ifdef MEMDUMMY

void get_procfs(p_vm,p_rss)
float *p_vm, *p_rss;
{
	*p_vm=0.0; *p_rss=0.0;
}

#else

/* get the values and return them */
void get_procfs(p_vm,p_rss)
float *p_vm, *p_rss;
{
    int ret = 0;
    float vsize,rss;
    float vsizemb, rssmb;
    ps_procstat_t p;

    if ( ( ret = ps_procstat(&p) ) == 0 ) {
	vsize = p.vsize;
	rss = p.rss;
    }

    pgsize = getpagesize();
    vsizemb = ( vsize ) / (float) MEGABYTE;
    rssmb = ( rss * pgsize ) / (float) MEGABYTE;

    *p_vm = vsizemb;
    *p_rss = rssmb;
}

#endif

int ps_procstat(ps_procstat_t *p)
{
    FILE *fd;
    char fname[64];

    sprintf(fname, "/proc/self/stat");

    if ( ( fd=fopen(fname, "r") ) == 0 ) {
        return -1;
    }

    fscanf(fd, "%d (%[^)]) %c %d", &p->pid, p->comm, &p->state, &p->ppid);
    fscanf(fd, "%d %d %d %d %lu",  &p->pgrp, &p->session, &p->tty, &p->tpgid,
           &p->flags);
    fscanf(fd, "%lu %lu %lu %lu",  &p->minflt, &p->cminflt, &p->majflt,
           &p->cmajflt);
    fscanf(fd, "%ld %ld %ld %ld",  &p->utime, &p->stime, &p->cutime,
           &p->cstime);
    fscanf(fd, "%ld %ld %*d %lu",  &p->priority, &p->nice, &p->itrealvalue);
    fscanf(fd, "%lu %lu %lu %lu",  &p->starttime, &p->vsize, &p->rss,
           &p->rlim);
    fscanf(fd, "%lu %lu %lu %lu %lu", &p->startcode, &p->endcode,
           &p->startstack, &p->kstkesp, &p->kstkeip);
    fscanf(fd, "%d %d %d %d %lu",   &p->signal, &p->blocked, &p->sigignore,
           &p->sigcatch, &p->wchan);
    fscanf(fd, "%lu %lu %d %d",    &p->nswap, &p->cnswap, &p->exitsignal,
           &p->processor);

    fclose(fd);

    return 0;
}

