#ifndef JABS_DEBUG_H
#define JABS_DEBUG_H
#include <stdio.h>

#ifdef DEBUG
    #ifdef __FILE_NAME__
    #define DEBUGMSG(format, ...) fprintf(stderr, "DEBUG %s [%s:%i] " format "\n", __func__, __FILE_NAME__, __LINE__, __VA_ARGS__);
    #define DEBUGSTR(str) fprintf(stderr, "DEBUG %s [%s:%i] %s\n", __func__, __FILE_NAME__, __LINE__, str);
    #else
    #define DEBUGMSG(format, ...) fprintf(stderr, "DEBUG %s [%s:%i] " format "\n", __func__, __FILE__, __LINE__, __VA_ARGS__);
    #define DEBUGSTR(str) fprintf(stderr, "DEBUG %s [%s:%i] %s\n", __func__, __FILE__, __LINE__, str);
#endif
#else
    #define DEBUGMSG(x, ...)
    #define DEBUGSTR(x)
#endif

#ifdef DEBUG_VERBOSE
#define DEBUGVERBOSEMSG(format, ...) DEBUGMSG(format, __VA_ARGS__)
#define DEBUGVERBOSESTR(str) DEBUGSTR(str)
#else
#define DEBUGVERBOSEMSG(format, ...)
#define DEBUGVERBOSESTR(str)
#endif

#endif // JABS_DEBUG_H
