#ifndef PROFILING_JSON_H_
#define PROFILING_JSON_H_

#define B4S_TRACING 1

#include <stdio.h>
#include <stdlib.h>

#define UNUSED(x) (void)(x)

#define STRINGIZE_DETAIL(x) #x
#define STRINGIZE(x) STRINGIZE_DETAIL(x)

#define EXPECT(cond, msg) \
    do { \
        if (!(cond)) { \
            fprintf(stderr, "FATAL ERROR: %s\n", msg); \
            exit(EXIT_FAILURE); \
        } \
    } while (0)

void profiling_json_init(void);
void profiling_json_finalize(void);
void profiling_json_emit(const char *cat, const char *name, const char *ph);

#ifdef B4S_TRACING

#include <stdint.h>

#define B4S_TRACE_BEGIN(cat, name) \
    profiling_json_emit(STRINGIZE(cat), name, "B")

#define B4S_TRACE_END(cat) \
    profiling_json_emit(#cat, "", "E")

#define B4S_TRACE_INSTANT(cat, name) \
    profiling_json_emit(#cat, name, "i")

#else // !B4S_TRACING

#define B4S_TRACE_BEGIN(cat, name) do { UNUSED(#cat); UNUSED(name); } while(0)
#define B4S_TRACE_END(cat) do { UNUSED(#cat); } while(0)
#define B4S_TRACE_INSTANT(cat, name) do { UNUSED(#cat); UNUSED(name); } while(0)

#endif // B4S_TRACING

#endif // PROFILING_JSON_H_
