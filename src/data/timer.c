/**
 * @file timer.c
 *
 * @brief maintain timers
 * copied from bamps 3.0 - slightly modified
*/

//#include "bamps.h"
#include "ader.h"

#include <sys/resource.h>
#include <sys/time.h>


/**
 * @struct tTimer
 * @brief struct represents a timer in bamps
 */
typedef struct {
  char *name;                   ///< Name of the timer
  int l;                        ///< Needed for sorting timer. Seems always to be -1
  double start;                 ///< Starting time
  double stop;                  ///< Ending time
  double time;                  ///< Difference between stop and start
  int n;                        ///< Counts how often counter has been started and stopped
} tTimer;

tTimer *tdb;                    ///< Timer database
int ntimers = 0;                ///< Number of active timers

static int timer_on = 0;        ///< Flag stating if timing is on or not




/**
 * @brief enable/disable timer
 */
void timer_enable(int on)
{
  timer_on = on;
}




/**
 * @brief return user time + system time
 *
 * add them up since some gpu calls show up as system time
 */
double time_user()
{
  struct rusage usage;
  long int sec, usec;
  double t;

  getrusage(RUSAGE_SELF, &usage);
  sec = usage.ru_utime.tv_sec + usage.ru_stime.tv_sec;
  usec = usage.ru_utime.tv_usec + usage.ru_stime.tv_usec;
  t = sec + usec * 1e-6;
  if (0)
    printf("%f\n", t);
  return t;
}




/**
 * @brief return wall time in seconds
 */
double time_wall()
{
#ifdef MPI
  return MPI_Wtime();
#else
  struct timeval tv;
  gettimeofday(&tv, 0);
  return tv.tv_sec + tv.tv_usec * 1e-6;
#endif
}




/**
 * @brief return time since last call
 */
double delta_time(void)
{
  static double t0 = -1;
  double dt;

  if (t0 < 0)
    t0 = time_wall();
  dt = fabs(time_wall() - t0);
  t0 += dt;
  return dt;
}




double timer[100];

/**
 * @brief add time in seconds since last call
 */
void add_time(int n)
{
  static double t0 = -1;
  double t;
  int i;

  if (n < 0 || n >= 100) {
    t0 = time_wall();
    for (i = 0; i < 100; i++)
      timer[i] = 0;
  } else {
    t = time_wall();
    timer[n] += t - t0;
    t0 = t;
  }
}



/**********************************************************************/


/**
 * @brief return pointer to entry for name, if name does not exist yet, allocate
*/
tTimer *timer_get( char *name )
{
  int l = -1;
  int i;

  for (i = 0; i < ntimers; i++)
    if (l == tdb[i].l && !strcmp(tdb[i].name, name))
      return tdb + i;

  // quick fix to let this run with threads
  // tdb = realloc(tdb, (++ntimers) * sizeof(tTimer));
  if (ntimers == 0)
    tdb = calloc(100, sizeof(tTimer));
  if (ntimers >= 100){
    //errorexit("fixme");
  }
  ntimers++;
  // fixme: not good enough, memory ok, but initialization might happen
  //   simultaneously, should initialize all timers serially, or lock

  tdb[i].name = strdup(name);
  tdb[i].l = l;
  tdb[i].start = -1;
  tdb[i].time = 0;
  tdb[i].n = 0;
  return tdb + i;
}




/**
 * @brief return time stored in timer
 */
double timer_time(char *name)
{
  return timer_get(name)->time;
}




/**
 * @brief start timer
 */
int timer_start(char *name)
{
  if (timer_on) {
    tTimer *t;

    if (timer_on == -1) {
      if (!timer_on)
        return 0;
    }

    t = timer_get(name);
    t->start = time_wall();

    if (0)
      printf("starting timer %s, %.15e\n", t->name, t->start);
  }
  return 1;
}




/**
 * @brief stop timer
 */
int timer_stop(char *name)
{
  if (timer_on) {
    tTimer *t = timer_get(name);

    if (t->start < 0)
      return 0;
    t->n += 1;
    t->time += time_wall() - t->start;
    t->start = -1;

    if (0)
      printf("stopping timer %s, %e\n", t->name, t->time);
  }
  return 1;
}




/**
 * @brief timer comparison for sorting
 */
int timer_compare(const void *a, const void *b)
{
  tTimer *t = (tTimer *) a;
  tTimer *u = (tTimer *) b;
  int i;

  if (t->l == -1 && u->l >= 0)
    return 1;
  if (t->l >= 0 && u->l == -1)
    return -1;

  if (t->l == -1 && u->l == -1) {
    if (t->time < u->time)
      return -1;
    if (t->time > u->time)
      return 1;
    return 0;
  }

  i = strcmp(t->name, u->name);
  if (!i) {
    if (t->l < u->l)
      return -1;
    if (t->l > u->l)
      return 1;
  }

  return i;
}




/**
 * @brief print timer info
 */
void timer_print( PARA *par )
{
//  tG *g = global_grid[0];
  FILE *fp;
  char *outdir = "../data" ; //getString("output.dir");
  char f[100], s[1000];
  tTimer *t;
  double p, total;
  int i;

  if (!timer_on)
    return;

  /* open file */
  snprintf(f, 100, "%%s/timer.t");
  snprintf(s, 1000, f, outdir);
  fp = fopen(s, "a");
  if (!fp) {
    printf("could not open %s", s);
    return;
  }
//  fprintf(fp, "\nTimers after evolution step %d, time %.3f\n",
//          g->evolve_step, g->time);

  fprintf ( fp, "\n#################################################\nTimers\n" );
  fprintf ( fp, "equation:\t%s\n", par->project );
  fprintf ( fp, "predictor:\t%s\n", par->scheme );
  fprintf ( fp, "M=%d, N=%d\n\n", par->M, par->N );

  /* update timer for main() */
  timer_stop("main");
  timer_start("main");
  t = timer_get("main");
  total = t->time;
  t->n = 1;

  /* sort */
  if (1)
    qsort(tdb, ntimers, sizeof(tTimer), timer_compare);

  /* print */
  for (i = 0; i < ntimers; i++) {
    t = tdb + i;
    p = t->time / total * 100;

    if (t->l < 0)
      fprintf(fp, "%-24s %5.1f  %9.4f %7d\n", t->name, p, t->time, t->n);
  }

  /* done */
  fclose(fp);
}




/**
 * @brief return time since last call
 *
 *  independent of detailed timers above
 *  used for M per seconds in log.t \n
*/
double time_since_last()
{
  static double lasttime = 0.0;
  if (lasttime == 0)
    lasttime = time_wall();
  double newtime = time_wall();
  double dt = newtime - lasttime;
  lasttime = newtime;
  return fabs(dt);
}


