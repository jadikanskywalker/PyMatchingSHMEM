#include "profiling_json.h"
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#ifdef USE_SHMEM
#include <shmem.h>
#endif
#include <sys/time.h>
#include <time.h>

#ifdef B4S_TRACING

static FILE *trace_file = NULL;
static int my_pe = -1;

static uint64_t get_time_us() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (uint64_t)ts.tv_sec * 1000000 + ts.tv_nsec / 1000;
}

void profiling_json_init(void) {
    char filename[64];
#ifdef USE_SHMEM
    my_pe = shmem_my_pe();
    snprintf(filename, sizeof(filename), "trace_%d.json", my_pe);
#else
    snprintf(filename, sizeof(filename), "trace.json");
#endif
    trace_file = fopen(filename, "w");
    EXPECT(trace_file, "failed opening trace file");

    // make sure trace written even if crash
    setvbuf(trace_file, NULL, _IOLBF, 0);

    fprintf(trace_file, "[\n");
    
    fprintf(trace_file, 
        "{\"name\":\"process_name\",\"ph\":\"M\",\"pid\":%d,\"args\":{\"name\":\"PE %d\"}},\n",
        my_pe, my_pe);
}

void profiling_json_finalize(void) {
    if (trace_file) {
        // dummy event closes json validly
        fprintf(trace_file, "{\"name\":\"root\",\"ph\":\"i\",\"pid\":%d,\"tid\":0,\"ts\":%llu}\n]\n", 
                my_pe, (unsigned long long)get_time_us());
        fclose(trace_file);
        trace_file = NULL;
    }
}

void profiling_json_emit(const char *cat, const char *name, const char *ph) {
    // SHOULD never happen, but i want to keep this non panicking
    if (!trace_file) return;

    uint64_t now = get_time_us();
    
    fprintf(trace_file, 
        "{\"cat\":\"%s\",\"name\":\"%s\",\"ph\":\"%s\",\"pid\":%d,\"tid\":%d,\"ts\":%llu},\n",
        cat, name, ph, my_pe, omp_get_thread_num(), (unsigned long long)now);
}

#else

void profiling_json_init(void) {}
void profiling_json_finalize(void) {}
void profiling_json_emit(const char *cat, const char *name, const char *ph) {
    UNUSED(cat);
    UNUSED(name);
    UNUSED(ph);
}

#endif
