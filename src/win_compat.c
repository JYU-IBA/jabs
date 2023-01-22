#include "win_compat.h"
#ifdef WIN32
#include <string.h>
char *strsep(char **stringp, const char *delim) {
	char *start= *stringp;
	char *p;

	p = (start != NULL) ? strpbrk(start, delim) : NULL;

	if(p==NULL) {
		*stringp=NULL;
	} else {
		*p = '\0';
		*stringp=p+1;
	}
	return start;
}

char *strndup(char *s1, size_t n) {
    char *buf = (char *) malloc(n + 1);
    if(!buf) {
        return 0;
    }
    size_t i;
    for(i = 0; ((i < n) && (s1[i] != 0)); i++) {
        buf[i] = s1[i];
    }
    buf[i] = 0;
    return buf;
}

#include <stdlib.h>
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <pathcch.h>
char *dirname(char *path) {
    size_t s = strlen(path)+1;
    size_t convertedChars = 0;
    wchar_t *pszPath = malloc(sizeof(wchar_t)*s);
    mbstowcs_s(&convertedChars, pszPath, s, path, _TRUNCATE); /* convert multibyte to wide */
    if(!pszPath)
        return NULL;
    HRESULT result = PathCchRemoveFileSpec(pszPath, convertedChars); /* remove trailing file name or directory */
    if(result != S_OK)
        return NULL;
    size_t s_out = 2*(wcslen(pszPath)+1);
    char *path_out = malloc(s_out);
    wcstombs_s(&convertedChars, path_out, s_out, pszPath, _TRUNCATE);
    free(pszPath);
    return path_out; /* Could be NULL */
}
#endif

#ifdef _MSC_VER
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>

int vasprintf(char **strp, const char *format, va_list ap)
{
    int len = vscprintf(format, ap);
    if (len == -1)
        return -1;
    char *str = (char*)malloc((size_t) len + 1);
    if (!str)
        return -1;
    int retval = vsnprintf(str, len + 1, format, ap);
    if (retval == -1) {
        free(str);
        return -1;
    }
    *strp = str;
    return retval;
}

int asprintf(char **strp, const char *format, ...)
{
    va_list ap;
    va_start(ap, format);
    int retval = vasprintf(strp, format, ap);
    va_end(ap);
    return retval;
}


#include <BaseTsd.h>
#define ssize_t SSIZE_T

/*	$NetBSD: getdelim.c,v 1.2 2015/12/25 20:12:46 joerg Exp $	*/
/*	NetBSD-src: getline.c,v 1.2 2014/09/16 17:23:50 christos Exp 	*/

/*-
 * Copyright (c) 2011 The NetBSD Foundation, Inc.
 * All rights reserved.
 *
 * This code is derived from software contributed to The NetBSD Foundation
 * by Christos Zoulas.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE NETBSD FOUNDATION, INC. AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE FOUNDATION OR CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include <stdio.h>
#include <stdlib.h>

ssize_t getdelim(char **buf, size_t *bufsiz, int delimiter, FILE *fp) {
    char *ptr, *eptr;
    if (*buf == NULL || *bufsiz == 0) {
        *bufsiz = BUFSIZ;
        if ((*buf = malloc(*bufsiz)) == NULL)
            return -1;
    }

    for (ptr = *buf, eptr = *buf + *bufsiz;;) {
        int c = fgetc(fp);
        if (c == -1) {
            if (feof(fp)) {
                ssize_t diff = (ssize_t) (ptr - *buf);
                if (diff != 0) {
                    *ptr = '\0';
                    return diff;
                }
            }
            return -1;
        }
        *ptr++ = c;
        if (c == delimiter) {
            *ptr = '\0';
            return ptr - *buf;
        }
        if (ptr + 2 >= eptr) {
            char *nbuf;
            size_t nbufsiz = *bufsiz * 2;
            ssize_t d = ptr - *buf;
            if ((nbuf = realloc(*buf, nbufsiz)) == NULL)
                return -1;
            *buf = nbuf;
            *bufsiz = nbufsiz;
            eptr = nbuf + nbufsiz;
            ptr = nbuf + d;
        }
    }
}

/*	$NetBSD: getline.c,v 1.2 2015/12/25 20:12:46 joerg Exp $	*/
/*	NetBSD-src: getline.c,v 1.2 2014/09/16 17:23:50 christos Exp	*/

ssize_t getline(char **buf, size_t *bufsiz, FILE *fp) {
    return getdelim(buf, bufsiz, '\n', fp);
}

#endif
