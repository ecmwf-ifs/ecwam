#include <stdint.h>
#include <stddef.h>

/* Byte-reverse every element of buf in-place, for converting restart
 * payload between native and big-endian order in one pass.
 * bswap32_array is for 4-byte reals, bswap64_array for 8-byte ones;
 * reversing a double is not the same as reversing its two halves. */

void bswap32_array(float *buf, int n)
{
    uint32_t *p = (uint32_t *)buf;
    size_t    i, len = (size_t)n;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (i = 0; i < len; i++)
        p[i] = __builtin_bswap32(p[i]);
}

void bswap64_array(double *buf, int n)
{
    uint64_t *p = (uint64_t *)buf;
    size_t    i, len = (size_t)n;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (i = 0; i < len; i++)
        p[i] = __builtin_bswap64(p[i]);
}
