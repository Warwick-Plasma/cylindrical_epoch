/*
 * SDF (Self-Describing Format) Software Library
 * Copyright (c) 2013-2016, SDF Development Team
 *
 * Distributed under the terms of the BSD 3-clause License.
 * See the LICENSE file for details.
 */

#ifndef _SDF_VECTOR_TYPE_H_
#define _SDF_VECTOR_TYPE_H_
#include <stdlib.h>
#include <string.h>

/* Growable array datatype and helper routines */

typedef struct vector_type vector_t;

struct vector_type {
    int *data;
    int allocated, size;
};

static vector_t *vector_new(void)
{
    vector_t *vector;

    vector = (vector_t*)malloc(sizeof(vector_t));
    vector->allocated = 32;
    vector->size = 0;
    vector->data = (int*)malloc(vector->allocated * sizeof(*vector->data));

    return vector;
}

static void vector_push_back(vector_t *vector, int val)
{
    int *data;

    // Grow vector if necessary
    if (vector->size == vector->allocated) {
        vector->allocated = vector->allocated << 1;
        data = (int*)malloc(vector->allocated * sizeof(*data));
        memcpy(data, vector->data, vector->size * sizeof(*data));
        free(vector->data);
        vector->data = data;
    }

    vector->data[vector->size++] = val;
}

static void vector_truncate(vector_t *vector)
{
    int *data;

    vector->allocated = vector->size;
    data = (int*)malloc(vector->allocated * sizeof(*data));
    memcpy(data, vector->data, vector->size * sizeof(*data));
    free(vector->data);
    vector->data = data;
}

static void vector_free(vector_t *vector)
{
    free(vector->data);
    free(vector);
}

#endif
