/*
 * SDF (Self-Describing Format) Software Library
 * Copyright (c) 2010-2016, SDF Development Team
 *
 * Distributed under the terms of the BSD 3-clause License.
 * See the LICENSE file for details.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sdf.h>
#include "sdf_input.h"
#include "sdf_input_cartesian.h"
#include "sdf_control.h"

#ifndef PARALLEL
# ifdef _WIN32
#  include <io.h>
#  define FSEEKO _fseeki64
# else
#  include <unistd.h>
#  include <sys/mman.h>
#  define FSEEKO fseeko
# endif
#endif

//#define SDF_COMMON_MESH_LENGTH (4 + 8 + h->id_length + 4 * b->ndims)

#define SDF_COMMON_MESH_INFO() do { \
    if (!h->current_block || !h->current_block->done_header) { \
        if (h->rank == h->rank_master) { \
            fprintf(stderr, "*** ERROR ***\n"); \
            fprintf(stderr, "SDF block header has not been read." \
                    " Ignoring call.\n"); \
        } \
        return 1; \
    } \
    b = h->current_block; \
    if (b->done_info) return 0; \
    h->current_location = b->block_start + h->block_header_length; \
    b->done_info = 1; } while(0)


int sdf_read_plain_mesh_info(sdf_file_t *h)
{
    sdf_block_t *b;
    int i;
    int32_t dims_in[SDF_MAXDIMS];
    int32_t *dims_ptr = dims_in;

    // Metadata is
    // - mults     REAL(r8), DIMENSION(ndims)
    // - labels    CHARACTER(id_length), DIMENSION(ndims)
    // - units     CHARACTER(id_length), DIMENSION(ndims)
    // - geometry  INTEGER(i4)
    // - minval    REAL(r8), DIMENSION(ndims)
    // - maxval    REAL(r8), DIMENSION(ndims)
    // - dims      INTEGER(i4), DIMENSION(ndims)

    SDF_COMMON_MESH_INFO();

    SDF_READ_ENTRY_ARRAY_REAL8(b->dim_mults, b->ndims);

    SDF_READ_ENTRY_ARRAY_ID(b->dim_labels, b->ndims);

    SDF_READ_ENTRY_ARRAY_ID(b->dim_units, b->ndims);

    SDF_READ_ENTRY_INT4(b->geometry);

    SDF_READ_ENTRY_ARRAY_REAL8(b->extents, 2*b->ndims);

    SDF_READ_ENTRY_ARRAY_INT4(dims_ptr, b->ndims);
    b->nelements = 0;
    for (i = 0; i < b->ndims; i++) {
        b->dims[i] = dims_in[i];
        b->nelements += b->dims[i];
    }

    b->stagger = SDF_STAGGER_VERTEX;
    for (i = 0; i < b->ndims; i++) b->const_value[i] = 1;

    b->ndim_labels = b->ndim_units = b->ndims;

    // Calculate per block parallel factorisation
    // This will be fixed up later once we have the whole block list.
    sdf_factor(h);

    if (b->blocktype == SDF_BLOCKTYPE_PLAIN_MESH) {
        b->nelements = b->nelements_local = 0;
        for (i = 0; i < b->ndims; i++) {
            b->nelements += b->dims[i];
            b->nelements_local += b->local_dims[i];
        }
    } else {
        b->nelements = b->nelements_local = 1;
        for (i = 0; i < b->ndims; i++) {
            b->nelements *= b->dims[i];
            b->nelements_local *= b->local_dims[i];
        }
    }

    return 0;
}



int sdf_read_plain_variable_info(sdf_file_t *h)
{
    sdf_block_t *b;
    int i;
    int32_t dims_in[SDF_MAXDIMS];
    int32_t *dims_ptr = dims_in;

    // Metadata is
    // - mult      REAL(r8)
    // - units     CHARACTER(id_length)
    // - meshid    CHARACTER(id_length)
    // - dims      INTEGER(i4), DIMENSION(ndims)
    // - stagger   INTEGER(i4)

    SDF_COMMON_MESH_INFO();

    SDF_READ_ENTRY_REAL8(b->mult);

    SDF_READ_ENTRY_ID(b->units);

    SDF_READ_ENTRY_ID(b->mesh_id);

    SDF_READ_ENTRY_ARRAY_INT4(dims_ptr, b->ndims);
    b->nelements = 1;
    for (i = 0; i < b->ndims; i++) {
        b->dims[i] = dims_in[i];
        b->nelements *= b->dims[i];
    }

    SDF_READ_ENTRY_INT4(b->stagger);
    for (i = 0; i < b->ndims; i++) b->const_value[i] = (b->stagger & 1<<i);

    // Calculate per block parallel factorisation
    // This will be fixed up later once we have the whole block list.
    sdf_factor(h);

    return 0;
}



static int sdf_create_1d_distribution(sdf_file_t *h, int global, int local,
        int start)
{
#ifdef PARALLEL
    sdf_block_t *b = h->current_block;
    int lengths[3];
    MPI_Aint disp[3];
    MPI_Datatype types[3];

    lengths[0] = 1;
    lengths[1] = local;
    lengths[2] = 1;
    disp[0] = 0;
    disp[1] = start * SDF_TYPE_SIZES[b->datatype];
    disp[2] = global * SDF_TYPE_SIZES[b->datatype];
    types[0] = MPI_LB;
    types[1] = b->mpitype;
    types[2] = MPI_UB;

    MPI_Type_create_struct(3, lengths, disp, types, &b->distribution);
    MPI_Type_commit(&b->distribution);
#endif
    return 0;
}



static int sdf_plain_mesh_distribution(sdf_file_t *h)
{
#ifdef PARALLEL
    sdf_block_t *b = h->current_block;
    int n;
    int sizes[SDF_MAXDIMS], subsizes[SDF_MAXDIMS];

    for (n=0; n < b->ndims; n++) {
        b->dims[n] -= 2 * b->ng;
        b->local_dims[n] -= 2 * b->ng;
        sizes[n] = (int)b->dims[n];
        subsizes[n] = (int)b->local_dims[n];
    }

    // Get starts for creating subarray
    sdf_factor(h);

    MPI_Type_create_subarray(b->ndims, sizes, subsizes, b->starts,
            MPI_ORDER_FORTRAN, b->mpitype, &b->distribution);
    MPI_Type_commit(&b->distribution);

    b->nelements_local = 1;
    for (n=0; n < b->ndims; n++) {
        b->dims[n] += 2 * b->ng;
        b->local_dims[n] += 2 * b->ng;
        b->nelements_local *= b->local_dims[n];
    }
#endif

    return 0;
}



static int sdf_free_distribution(sdf_file_t *h)
{
#ifdef PARALLEL
    sdf_block_t *b = h->current_block;

    if (b->ng) return 0;

    MPI_Type_free(&b->distribution);
#endif
    return 0;
}



static int sdf_helper_read_array_halo(sdf_file_t *h, void **var_in)
{
    sdf_block_t *b = h->current_block;
    char **var_ptr = (char**)var_in;
    char *read_var = *var_ptr;
    char convert;
    int sz;
    size_t count;
#ifdef PARALLEL
    MPI_Datatype distribution, facetype;
    int sizes[SDF_MAXDIMS], subsizes[SDF_MAXDIMS];
    int face[SDF_MAXDIMS];
    int i, tag;
    int64_t offset;
    char *p1, *p2;
#else
    char *read_ptr;
    int i, j, k;
    int nx, ny, nz;
    int j0, j1, k0, k1;
    int subsize;
#endif

    count = b->nelements_local;

    if (*var_ptr) free(*var_ptr);
    sz = SDF_TYPE_SIZES[b->datatype];
    *var_ptr = read_var = calloc(1, count * sz);

#ifdef PARALLEL
    for (i=0; i < b->ndims; i++) {
        sizes[i] = b->local_dims[i];
        subsizes[i] = sizes[i] - 2 * b->ng;
        b->starts[i] = b->ng;
    }
    for (i=b->ndims; i < SDF_MAXDIMS; i++)
        subsizes[i] = 1;

    MPI_Type_create_subarray(b->ndims, sizes, subsizes, b->starts,
            MPI_ORDER_FORTRAN, b->mpitype, &distribution);

    MPI_Type_commit(&distribution);

    MPI_File_set_view(h->filehandle, h->current_location, b->mpitype,
            b->distribution, "native", MPI_INFO_NULL);
    MPI_File_read_all(h->filehandle, read_var, 1, distribution,
            MPI_STATUS_IGNORE);
    MPI_File_set_view(h->filehandle, 0, MPI_BYTE, MPI_BYTE, "native",
            MPI_INFO_NULL);

    MPI_Type_free(&distribution);

    // Swap ghostcell faces
    for (i=0; i < b->ndims; i++) {
        face[i] = sizes[i] - 2 * b->ng;
        b->starts[i] = b->ng;
    }

    tag = 1;
    offset = sz;
    for (i=0; i < b->ndims; i++) {
        face[i] = b->ng;
        b->starts[i] = 0;

        MPI_Type_create_subarray(b->ndims, sizes, face, b->starts,
                MPI_ORDER_FORTRAN, b->mpitype, &facetype);
        MPI_Type_commit(&facetype);

        p1 = (char*)b->data + b->ng * offset;
        p2 = (char*)b->data + (sizes[i] - b->ng) * offset;
        MPI_Sendrecv(p1, 1, facetype, b->proc_min[i], tag, p2, 1, facetype,
                b->proc_max[i], tag, h->comm, MPI_STATUS_IGNORE);
        tag++;

        p1 = (char*)b->data + (sizes[i] - 2 * b->ng) * offset;
        p2 = b->data;
        MPI_Sendrecv(p1, 1, facetype, b->proc_max[i], tag, p2, 1, facetype,
                b->proc_min[i], tag, h->comm, MPI_STATUS_IGNORE);
        tag++;

        MPI_Type_free(&facetype);

        face[i] = sizes[i];
        offset *= b->local_dims[i];
    }
#else
    j0 = k0 = 0;
    j1 = k1 = ny = nz = 1;

    i = b->ng;
    nx = b->local_dims[0];
    subsize = nx - 2 * i;

    if (b->ndims > 1) {
        ny = b->local_dims[1];
        j0 = b->ng;
        j1 = ny - b->ng;
    }
    if (b->ndims > 2) {
        nz = b->local_dims[2];
        k0 = b->ng;
        k1 = nz - b->ng;
    }

    for (k = k0; k < k1; k++) {
    for (j = j0; j < j1; j++) {
        read_ptr = read_var + (i + nx * (j + ny * k)) * sz;
        FSEEKO(h->filehandle, h->current_location, SEEK_SET);
        if (!fread(read_ptr, sz, subsize, h->filehandle)) return 1;
        h->current_location += subsize * sz;
    }}
#endif

    if (h->use_float && b->datatype == SDF_DATATYPE_REAL8)
        convert = 1;
    else
        convert = 0;

    if (convert) {
        size_t i;
        float *r4;
        double *old_var, *r8;
        r8 = old_var = (double*)read_var;
        r4 = (float*)(*var_ptr);
        if (!r4) {
            *var_ptr = malloc(count * sizeof(float));
            r4 = (float*)(*var_ptr);
        }
        for (i=0; i < count; i++)
            *r4++ = (float)(*r8++);
        if (!h->mmap) free(old_var);
        b->datatype_out = SDF_DATATYPE_REAL4;
#ifdef PARALLEL
        b->mpitype_out = MPI_FLOAT;
#endif
    }

    return (count * SDF_TYPE_SIZES[b->datatype_out]);
}



int64_t sdf_helper_read_array(sdf_file_t *h, void **var_in, int dim)
{
    sdf_block_t *b = h->current_block;
    char **var_ptr = (char**)var_in;
    char *read_ptr = *var_ptr, *read_var = *var_ptr;
    char convert;
    int i, sz;
    size_t count, length;
#ifdef PARALLEL
    int64_t dims[SDF_MAXDIMS] = {0};
    int64_t starts[SDF_MAXDIMS] = {0};
    int64_t ends[SDF_MAXDIMS] = {1,1,1,1};
    size_t nelements;
#else
    size_t mlen, mstart, moff;
    int64_t *offset_starts = NULL, *offset_ends = NULL;
    int64_t ncount, offset, nreads;
    int *loop_counts = NULL, *idx = NULL;
    int n, ndims, combined_reads = 0;
#endif

    if (b->ndims < 1) return 0;

    if (b->ng) return sdf_helper_read_array_halo(h, var_in);

    if (dim < 0) {
        count = 1;
        for (i = 0; i < b->ndims; i++)
            count *= b->local_dims[i];
    } else
        count = b->local_dims[dim];

#ifdef PARALLEL
    nelements = count;
#endif

    if (b->array_starts) {
        if (dim < 0) {
            count = 1;
            for (i = 0; i < b->ndims; i++)
                count *= (b->array_ends[i] - b->array_starts[i]);
        } else
            count = b->array_ends[dim] - b->array_starts[dim];
    }

    sz = SDF_TYPE_SIZES[b->datatype];
    length = sz * count;

#if !defined(PARALLEL) && !defined(_WIN32)
    if (h->mmap) {
        mlen = sysconf(_SC_PAGESIZE);
        mstart = mlen * (h->current_location / mlen);
        moff = h->current_location - mstart;
        b->mmap_len = mlen = length + moff;
        b->mmap = mmap(NULL, mlen, PROT_READ, MAP_SHARED, h->fd, mstart);
        if (*var_ptr) free(*var_ptr);
        *var_ptr = moff + b->mmap;
        return length;
    }
#endif

    if (h->use_float && b->datatype == SDF_DATATYPE_REAL8) {
        convert = 1;
        read_ptr = read_var = malloc(length);
    } else {
        convert = 0;
        if (!read_var) *var_ptr = read_ptr = read_var = malloc(length);
    }

#ifdef PARALLEL
    if (count != nelements)
        read_ptr = malloc(sz * nelements);

    MPI_File_set_view(h->filehandle, h->current_location, b->mpitype,
            b->distribution, "native", MPI_INFO_NULL);
    MPI_File_read_all(h->filehandle, read_ptr, nelements, b->mpitype,
            MPI_STATUS_IGNORE);
    MPI_File_set_view(h->filehandle, 0, MPI_BYTE, MPI_BYTE, "native",
            MPI_INFO_NULL);

    if (count != nelements && dim < 0) {
        memcpy(dims, b->local_dims, b->ndims * sizeof(*dims));
        memcpy(starts, b->array_starts, b->ndims * sizeof(*starts));
        memcpy(ends, b->array_ends, b->ndims * sizeof(*ends));

        if (b->datatype == SDF_DATATYPE_INTEGER4
                || b->datatype == SDF_DATATYPE_REAL4) {
            size_t i, j, k, l;
            int32_t *v1 = (int32_t*)*var_ptr;
            int32_t *v2 = (int32_t*)read_ptr;
            for (l = starts[3]; l < ends[3]; l++) {
            for (k = starts[2]; k < ends[2]; k++) {
            for (j = starts[1]; j < ends[1]; j++) {
            for (i = starts[0]; i < ends[0]; i++) {
                *v1 = v2[i + dims[0] * (j + dims[1] * (k + dims[2] * l))];
                v1++;
            }}}}
        } else if (b->datatype == SDF_DATATYPE_INTEGER8
                || b->datatype == SDF_DATATYPE_REAL8) {
            size_t i, j, k, l;
            int64_t *v1 = (int64_t*)*var_ptr;
            int64_t *v2 = (int64_t*)read_ptr;
            for (l = starts[3]; l < ends[3]; l++) {
            for (k = starts[2]; k < ends[2]; k++) {
            for (j = starts[1]; j < ends[1]; j++) {
            for (i = starts[0]; i < ends[0]; i++) {
                *v1 = v2[i + dims[0] * (j + dims[1] * (k + dims[2] * l))];
                v1++;
            }}}}
        }
        free(read_ptr);
    } else if (count != nelements) {
        memcpy(dims, b->local_dims, b->ndims * sizeof(*dims));
        memcpy(starts, b->array_starts, b->ndims * sizeof(*starts));
        memcpy(ends, b->array_ends, b->ndims * sizeof(*ends));

        if (b->datatype == SDF_DATATYPE_INTEGER4
                || b->datatype == SDF_DATATYPE_REAL4) {
            size_t i;
            int32_t *v1 = (int32_t*)*var_ptr;
            int32_t *v2 = (int32_t*)read_ptr;
            for (i = b->array_starts[dim]; i < b->array_ends[dim]; i++) {
                *v1 = v2[i];
                v1++;
            }
        } else if (b->datatype == SDF_DATATYPE_INTEGER8
                || b->datatype == SDF_DATATYPE_REAL8) {
            size_t i;
            int64_t *v1 = (int64_t*)*var_ptr;
            int64_t *v2 = (int64_t*)read_ptr;
            for (i = b->array_starts[dim]; i < b->array_ends[dim]; i++) {
                *v1 = v2[i];
                v1++;
            }
        }
        free(read_ptr);
    }
#else
    nreads = 1;
    ndims = 0;

    if (b->array_starts && dim < 0) {
        // First check for any reads which can be combined
        length = sz;
        for (i = 0; i < b->ndims; i++) {
            ncount = b->array_ends[i] - b->array_starts[i];
            length *= ncount;
            if (ncount == b->dims[i]) {
                combined_reads++;
                sz *= b->dims[i];
            } else
                break;
        }

        nreads = 1;
        ndims = b->ndims - combined_reads;
        if (ndims > 0) {
            offset_starts = malloc(ndims * sizeof(*offset_starts));
            offset_ends   = malloc(ndims * sizeof(*offset_ends));
            loop_counts   = malloc(ndims * sizeof(*loop_counts));
            idx           = malloc(ndims * sizeof(*idx));

            n = combined_reads;
            for (i = 0; i < ndims; i++) {
                offset_starts[i] = b->array_starts[n] * sz;
                offset_ends[i] = (b->dims[n] - b->array_ends[n]) * sz;
                sz *= b->dims[n];
                n++;
                if (n < ndims)
                    loop_counts[i] = b->array_ends[n] - b->array_starts[n];
                else
                    loop_counts[i] = 1;
                nreads *= loop_counts[i];
                idx[i] = 0;
            }
        }
    } else if (b->array_starts && dim >= 0) {
        ndims = 1;
        offset_starts = malloc(ndims * sizeof(*offset_starts));
        offset_ends   = malloc(ndims * sizeof(*offset_ends));
        loop_counts   = calloc(ndims, sizeof(*loop_counts));
        idx           = calloc(ndims, sizeof(*idx));
        offset_starts[0] = b->array_starts[dim] * sz;
        offset_ends[0] = (b->dims[dim] - b->array_ends[dim]) * sz;
        length = (b->array_ends[dim] - b->array_starts[dim]) * sz;
    }

    // Loop over an arbitrary number of dimensions reading array chunks

    offset = h->current_location;

    for (n = 0; n < nreads; n++) {
        for (i = 0; i < ndims; i++) {
            offset += offset_starts[i];
            if (idx[i] != 0)
                break;
        }

        FSEEKO(h->filehandle, offset, SEEK_SET);
        if (!fread(read_ptr, 1, length, h->filehandle))
            return (offset - h->current_location);

        read_ptr += length;
        offset += length;

        for (i = 0; i < ndims; i++) {
            offset += offset_ends[i];
            idx[i]++;
            if (idx[i] != loop_counts[i]) break;
            idx[i] = 0;
        }
    }

    if (ndims > 0) {
        free(offset_starts);
        free(offset_ends);
        free(loop_counts);
        free(idx);
    }
#endif
    if (h->swap) {
        if (b->datatype == SDF_DATATYPE_INTEGER4
                || b->datatype == SDF_DATATYPE_REAL4) {
            size_t i;
            int32_t *v = (int32_t*)*var_ptr;
            for (i=0; i < count; i++) {
                _SDF_BYTE_SWAP32(*v);
                v++;
            }
        } else if (b->datatype == SDF_DATATYPE_INTEGER8
                || b->datatype == SDF_DATATYPE_REAL8) {
            size_t i;
            int64_t *v = (int64_t*)*var_ptr;
            for (i=0; i < count; i++) {
                _SDF_BYTE_SWAP64(*v);
                v++;
            }
        }
    }

    if (convert) {
        size_t i;
        float *r4;
        double *old_var, *r8;
        r8 = old_var = (double*)read_var;
        r4 = (float*)(*var_ptr);
        if (!r4) {
            *var_ptr = malloc(count * sizeof(float));
            r4 = (float*)(*var_ptr);
        }
        for (i=0; i < count; i++)
            *r4++ = (float)(*r8++);
        if (!h->mmap) free(old_var);
        b->datatype_out = SDF_DATATYPE_REAL4;
#ifdef PARALLEL
        b->mpitype_out = MPI_FLOAT;
#endif
    }

    return (count * SDF_TYPE_SIZES[b->datatype_out]);
}



int sdf_read_plain_mesh(sdf_file_t *h)
{
    sdf_block_t *b = h->current_block;
    int n;

    if (b->done_data) return 0;
    if (!b->done_info) sdf_read_blocklist(h);

    h->current_location = b->data_location;

    if (!b->grids) {
        b->ngrids = 3;
        b->grids = calloc(b->ngrids, sizeof(float*));
    }

    if (h->print) {
        h->indent = 0;
        SDF_DPRNT("\n");
        SDF_DPRNT("b->name: %s ", b->name);
        for (n=0; n < b->ndims; n++) SDF_DPRNT("%" PRIi64 " ",b->local_dims[n]);
        SDF_DPRNT("\n");
        h->indent = 2;
    }

    for (n = 0; n < 3; n++) {
        if (b->ndims > n) {
#ifdef PARALLEL
            sdf_create_1d_distribution(h, (int)b->dims[n],
                    (int)b->local_dims[n], b->starts[n]);
#endif
            sdf_helper_read_array(h, &b->grids[n], n);
            sdf_free_distribution(h);
            if (h->print) {
                SDF_DPRNT("%s: ", b->dim_labels[n]);
                SDF_DPRNTar(b->grids[n], b->local_dims[n]);
            }
            h->current_location = h->current_location
                    + SDF_TYPE_SIZES[b->datatype] * b->dims[n];
        }
    }

    b->done_data = 1;

    return 0;
}



int sdf_read_lagran_mesh(sdf_file_t *h)
{
    sdf_block_t *b = h->current_block;
    int n;
    int64_t nelements = 1;

    if (b->done_data) return 0;
    if (!b->done_info) sdf_read_blocklist(h);

    h->current_location = b->data_location;

    if (!b->grids) {
        b->ngrids = 3;
        b->grids = calloc(b->ngrids, sizeof(float*));
    }

    if (h->print) {
        h->indent = 0;
        SDF_DPRNT("\n");
        SDF_DPRNT("b->name: %s ", b->name);
        for (n=0; n < b->ndims; n++) SDF_DPRNT("%" PRIi64 " ",b->local_dims[n]);
        SDF_DPRNT("\n");
        h->indent = 2;
    }

    sdf_plain_mesh_distribution(h);

    for (n = 0; n < b->ndims; n++) nelements *= b->dims[n];

    for (n = 0; n < 3; n++) {
        if (b->ndims > n) {
            sdf_helper_read_array(h, &b->grids[n], -1);
            sdf_convert_array_to_float(h, &b->grids[n], b->nelements_local);
            if (h->print) {
                SDF_DPRNT("%s: ", b->dim_labels[n]);
                SDF_DPRNTar(b->grids[n], b->nelements_local);
            }
            h->current_location = h->current_location
                    + SDF_TYPE_SIZES[b->datatype] * nelements;
        } else {
            b->grids[n] = calloc(1, SDF_TYPE_SIZES[b->datatype]);
        }
    }

    sdf_free_distribution(h);

    b->done_data = 1;

    return 0;
}



int sdf_read_plain_variable(sdf_file_t *h)
{
    sdf_block_t *b = h->current_block;
    int n;

    if (b->done_data) return 0;
    if (!b->done_info) sdf_read_plain_variable_info(h);

    h->current_location = b->data_location;

    sdf_plain_mesh_distribution(h);

    sdf_helper_read_array(h, &b->data, -1);

    sdf_free_distribution(h);

    if (h->print) {
        h->indent = 0;
        SDF_DPRNT("\n");
        SDF_DPRNT("b->name: %s ", b->name);
        for (n=0; n < b->ndims; n++) SDF_DPRNT("%" PRIi64 " ",b->local_dims[n]);
        SDF_DPRNT("\n  ");
        SDF_DPRNTar(b->data, b->nelements_local);
    }

    b->done_data = 1;

    return 0;
}
