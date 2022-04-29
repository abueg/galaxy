#pragma once

#include <stdio.h>

#ifdef TIMERS

void printTimers(FILE * out);

struct Timer;

#define TIMER_STRING(X) #X
#define TIMER_NAME(X) timer_ ## X ## _timer

#define START(X) startTimer(&(TIMER_NAME(X)))
#define STOP(X) stopTimer(&(TIMER_NAME(X)))

#define TIMER(X) extern struct Timer TIMER_NAME(X);
#include "timer_list.h"

void startTimer(struct Timer * timer);
void stopTimer(struct Timer * timer);

#else

static inline void printTimers(FILE * const file) {}

#define START(X)
#define STOP(X)
#define TIMER(X)

#endif
