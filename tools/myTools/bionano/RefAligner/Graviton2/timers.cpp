#ifdef TIMERS

#include <assert.h>
#include <omp.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>

#include "timers.h"

static int first = 1;

struct Timer
{
#ifdef DEBUG_TIMER_START_STOP
	char * name;
#endif
	pthread_mutex_t mutex;
	double elapsed;
	double before;
	double exclude;
	long nCalls;
	
	struct Timer * thread;
	struct Timer * caller;
	int nThreads;
	int callerLevel;
};

struct TimerStats
{
#ifdef DEBUG_PRINT_TIMERS
	const struct Timer * timer;
#endif
	char * name;
	double sumInclusive;
	double sumExclusive;
	double maxInclusive;
	double maxExclusive;
	double minInclusive;
	double minExclusive;
	long sumCalls;
	int nThreads;
};

struct TimerStack
{
	pthread_mutex_t mutex;
	struct Timer * current;
	
	struct TimerStack * parent;
	struct TimerStack * thread;
	int nThreads;
};

static struct TimerStack root;

#undef TIMER
#define TIMER(X) struct Timer TIMER_NAME(X);
#include "timer_list.h"

#ifdef DEBUG_TIMER_START_STOP
static void initTimer(struct Timer * const t, const char name[])
{
	t->name = strdup(name);
#else
static void initTimer(struct Timer * const t)
{
#endif
	pthread_mutex_init(&(t->mutex), NULL);
	t->elapsed = 0.0;
	t->before = 0.0;
	t->exclude = 0.0;
	t->nCalls = 0;
	t->thread = NULL;
	t->nThreads = 0;
	t->caller = NULL;
	t->callerLevel = -1;
}

static void initStack(struct TimerStack * const s, struct TimerStack * const parent)
{
	pthread_mutex_init(&(s->mutex), NULL);
	s->current = NULL;
	s->parent = parent;
	s->thread = NULL;
	s->nThreads = 0;
}

static void getStats(struct TimerStats * const s, const struct Timer * const t)
{
	if (t->nCalls) {
		assert(t->before == 0.0);
		const double inclusive = t->elapsed;
		const double exclusive = inclusive-(t->exclude);
		s->sumInclusive += inclusive;
		s->sumExclusive += exclusive;
		if (s->maxExclusive < exclusive) {
			s->maxExclusive = exclusive;
			s->maxInclusive = inclusive;
		}
		if (exclusive < s->minExclusive || s->minExclusive == 0.0) {
			s->minExclusive = exclusive;
			s->minInclusive = inclusive;
		}
		s->sumCalls += t->nCalls;
		++(s->nThreads);
	}

	for (int i = 0; i < t->nThreads; ++i) getStats(s, (t->thread)+i);
}

static void initStats(struct TimerStats * const s, const struct Timer * const t, const char name[])
{
#ifdef DEBUG_PRINT_TIMERS
	s->timer = t;
#endif
	s->name = strdup(name);
	s->sumInclusive = 0.0;
	s->sumExclusive = 0.0;
	s->maxInclusive = 0.0;
	s->maxExclusive = 0.0;
	s->minInclusive = 0.0;
	s->minExclusive = 0.0;
	s->sumCalls = 0;
	s->nThreads = 0;
	getStats(s, t);
};

static int compareStats(const void * const vleft, const void * const vright)
{
	const double left = ((const struct TimerStats *)vleft)->maxExclusive;
	const double right = ((const struct TimerStats *)vright)->maxExclusive;
	// Reverse order
	if (left < right) return 1;
	if (left > right) return -1;
	return 0;
}

static void printSeparator(FILE * const file, const int lineWidth)
{
	for (int i = 0; i < lineWidth; ++i) fputc('-', file);
	fputc('\n', file);
}

#ifdef DEBUG_PRINT_TIMERS
static void printTimerThreads(FILE * const file, const struct Timer * const t, const int level, const int threadNum[level])
{
	if (t->nThreads) {
		if (level > 0) {
			fprintf(file,"(%i", threadNum[0]);
			for (int i = 1; i < level; ++i) fprintf(file, ",%i", threadNum[i]);
			fprintf(file, ") ");
		}
		for (int i = 0; i < t->nThreads; ++i) {
			const double inclusive = t->thread[i].elapsed;
			const double exclusive = inclusive-(t->thread[i].exclude);
			fprintf(file, "%i", i);
			if (t->thread[i].nCalls) fprintf(file, "[%li,%lg,%lg]", t->thread[i].nCalls, exclusive, inclusive);
			fputc(' ', file);
		}
		fputc('\n', file);
		int newThreadNum[level+1];
		for (int j = 0; j < level; ++j) newThreadNum[j] = threadNum[j];
		for (int i = 0; i < t->nThreads; ++i) {
			newThreadNum[level] = i;
			printTimerThreads(file, (t->thread)+i, level+1, newThreadNum);
		}
	}
}
#endif

#undef TIMER
#ifdef DEBUG_TIMER_START_STOP
#define TIMER(X) initTimer(&(TIMER_NAME(X)), TIMER_STRING(X));
#else
#define TIMER(X) initTimer(&(TIMER_NAME(X)));
#endif

void startTimer(struct Timer * t)
{
	#pragma omp flush
	if (first) {
		#pragma omp flush
		#pragma omp critical(timers)
		if (first) {
			#include "timer_list.h"
			initStack(&root, NULL);
			#pragma omp flush
			first = 0;
			#pragma omp flush
		}
	}
	
	assert(t);
	const int level = omp_get_level();
	struct TimerStack * s = &root;

#ifdef DEBUG_TIMER_START_STOP
	#pragma omp critical
	{
		fprintf(stderr, "THREAD ");
		for (int l = 0; l <= level; ++l) fprintf(stderr, "%i ", omp_get_ancestor_thread_num(l));
		fprintf(stderr, " START %s\n", t->name);
	}
#endif
	
	// Find or create timer and stack for current thread
	for (int l = 1; l <= level; ++l) {
		const int atn = omp_get_ancestor_thread_num(l);
		#pragma omp flush
		if (!t->nThreads) {
			#pragma omp flush
			pthread_mutex_lock(&(t->mutex));
			if (!t->nThreads) {
				const int nThreads = omp_get_num_procs();
				t->thread = static_cast<Timer *>(malloc(nThreads*sizeof(*t)));
				assert(t->thread);
				for (int i = 0; i < nThreads; ++i) {
#ifdef DEBUG_TIMER_START_STOP
					initTimer((t->thread)+i, t->name);
#else
					initTimer((t->thread)+i);
#endif
				}
				#pragma omp flush
				t->nThreads = nThreads;
				#pragma omp flush
			}
			pthread_mutex_unlock(&(t->mutex));
		}
		assert(atn < t->nThreads);
		t = (t->thread)+atn;
		
		#pragma omp flush
		if (!s->nThreads) {
			#pragma omp flush
			pthread_mutex_lock(&(s->mutex));
			if (!s->nThreads) {
				const int nThreads = omp_get_num_procs();
				s->thread = static_cast<TimerStack *>(malloc(nThreads*sizeof(*s)));
				assert(s->thread);
				for (int i = 0; i < nThreads; ++i) initStack((s->thread)+i, s);
				#pragma omp flush
				s->nThreads = nThreads;
				#pragma omp flush
			}
			pthread_mutex_unlock(&(s->mutex));
		}
		s = (s->thread)+atn;
	}

	if (s->current) {
		t->caller = s->current;
		t->callerLevel = level;
		s->current = t;
	} else {
		s->current = t;
		t->caller = NULL;
		t->callerLevel = -1;
		for (int l = level; l > 0; --l) {
			if (omp_get_ancestor_thread_num(l) > 0) break;
			s = s->parent;
			assert(s);
			if (s->current) {
				t->caller = s->current;
				t->callerLevel = l-1;
				break;
			}
		}
	}
	
	++(t->nCalls);
	assert(t->before == 0.0);
	t->before = omp_get_wtime();
}

void stopTimer(struct Timer * t)
{
	const double after = omp_get_wtime();
	assert(!first);
	assert(t);

	const int level = omp_get_level();
	struct TimerStack * s = &root;
	
#ifdef DEBUG_TIMER_START_STOP
	#pragma omp critical
	{
		fprintf(stderr, "THREAD ");
		for (int l = 0; l <= level; ++l) fprintf(stderr, "%i ", omp_get_ancestor_thread_num(l));
		fprintf(stderr, " STOP %s\n", t->name);
	}
#endif
	
	// Find timer and stack for current thread
	for (int l = 1; l <= level; ++l) {
		const int atn = omp_get_ancestor_thread_num(l);
		assert(atn < t->nThreads);
		t = (t->thread)+atn;
		assert(atn < s->nThreads);
		s = (s->thread)+atn;
	}
	assert(t->before);
	const double elapsed = after-(t->before);
	t->elapsed += elapsed;
	assert(s->current == t);
	s->current = (t->callerLevel == level) ? t->caller : NULL;
	if (t->caller) t->caller->exclude += elapsed;
	t->before = 0.0;
	t->caller = NULL;
	t->callerLevel = -1;
}

void printTimers(FILE * const file)
{
	int nt = 0;
	#undef TIMER
	#define TIMER(X) ++nt;
	#include "timer_list.h"
	
	struct TimerStats s[nt];
	{
		int i = 0;
		#undef TIMER
		#define TIMER(X) initStats(s+i, &(TIMER_NAME(X)), TIMER_STRING(X)); ++i;
		#include "timer_list.h"
	}
	
	qsort(s, nt, sizeof(struct TimerStats), compareStats);
	
	int maxLen = strlen("Timer");
	for (int i = 0; i < nt; ++i) {
		const int len = s[i].sumCalls ? strlen(s[i].name) : 0;
		maxLen = (maxLen < len) ? len : maxLen;
	}

	const char title[] = "| %-*s | Thd | Call/T |  Max/T (Inclusive)  |  Avg/T (Inclusive)  |  Min/T (Inclusive)  |\n";
	const int lineWidth = maxLen+strlen(title)-5;
	printSeparator(file, lineWidth);
	fprintf(file, title, maxLen, "Timer");
	printSeparator(file, lineWidth);

	int unused = 0;
	for (int i = 0; i < nt; ++i) {
		if (s[i].sumCalls) {
			const double divThreads = 1.0/(s[i].nThreads);
			const double avgCalls = (s[i].sumCalls)*divThreads;
			const double avgExclusive = (s[i].sumExclusive)*divThreads;
			const double avgInclusive = (s[i].sumInclusive)*divThreads;
			fprintf(file, "| %-*s | %3i | %6.0lf | %8.2lf (%8.2lf) | %8.2lf (%8.2lf) | %8.2lf (%8.2lf) |\n", maxLen, s[i].name, s[i].nThreads, avgCalls, s[i].maxExclusive, s[i].maxInclusive, avgExclusive, avgInclusive, s[i].minExclusive, s[i].minInclusive);
		} else {
			++unused;
		}
	}
	printSeparator(file, lineWidth);
	if (unused) {
		fprintf(file, "Unused timers:");
		for (int i = 0; i < nt; ++i) if (s[i].sumCalls == 0) fprintf(file, " %s", s[i].name);
		fputc('\n', file);
	}
#ifdef DEBUG_PRINT_TIMERS
	printSeparator(file, 80);
	for (int i = 0; i < nt; ++i) {
		if (s[i].sumCalls) {
			fprintf(file, "%s ", s[i].name);
			const struct Timer * const t = s[i].timer;
			if (t->nCalls) {
				const double inclusive = t->elapsed;
				const double exclusive = inclusive-(t->exclude);
				fprintf(file, "[%li,%lg,%lg]", t->nCalls, exclusive, inclusive);
			}
			fputc('\n', file);
			printTimerThreads(file, t, 0, NULL);
			printSeparator(file, 80);
		}
	}
#endif
	fflush(file);
	for (int i = 0; i < nt; ++i) free(s[i].name);
}

#endif
