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
#include <sys/stat.h>
#include <assert.h>
#include <sdf.h>
#include "sdf_control.h"
#include "sdf_util.h"
#include "sdf_extension_util.h"
#include "commit_info.h"

#ifdef PARALLEL
# include <mpi.h>
#else
# ifndef _WIN32
#  include <sys/mman.h>
# endif
#endif

#ifdef _WIN32
# include <io.h>
# define FSEEKO _fseeki64
#else
# include <unistd.h>
# define FSEEKO fseeko
#endif

/**
 * @defgroup control
 * @brief Routines for controlling an SDF file
 */

#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define ABS(a) (((a) > 0) ? (a) : (-(a)))


#define SDF_MAX_RND (1LL<<32)

static uint32_t sdf_random(void);
static void sdf_random_init(void);
static int sdf_free_handle(sdf_file_t *h);
static int sdf_fclose(sdf_file_t *h);

const char *sdf_blocktype_c[] = {
    "SDF_BLOCKTYPE_NULL",
    "SDF_BLOCKTYPE_PLAIN_MESH",
    "SDF_BLOCKTYPE_POINT_MESH",
    "SDF_BLOCKTYPE_PLAIN_VARIABLE",
    "SDF_BLOCKTYPE_POINT_VARIABLE",
    "SDF_BLOCKTYPE_CONSTANT",
    "SDF_BLOCKTYPE_ARRAY",
    "SDF_BLOCKTYPE_RUN_INFO",
    "SDF_BLOCKTYPE_SOURCE",
    "SDF_BLOCKTYPE_STITCHED_TENSOR",
    "SDF_BLOCKTYPE_STITCHED_MATERIAL",
    "SDF_BLOCKTYPE_STITCHED_MATVAR",
    "SDF_BLOCKTYPE_STITCHED_SPECIES",
    "SDF_BLOCKTYPE_SPECIES",
    "SDF_BLOCKTYPE_PLAIN_DERIVED",
    "SDF_BLOCKTYPE_POINT_DERIVED",
    "SDF_BLOCKTYPE_CONTIGUOUS_TENSOR",
    "SDF_BLOCKTYPE_CONTIGUOUS_MATERIAL",
    "SDF_BLOCKTYPE_CONTIGUOUS_MATVAR",
    "SDF_BLOCKTYPE_CONTIGUOUS_SPECIES",
    "SDF_BLOCKTYPE_CPU_SPLIT",
    "SDF_BLOCKTYPE_STITCHED_OBSTACLE_GROUP",
    "SDF_BLOCKTYPE_UNSTRUCTURED_MESH",
    "SDF_BLOCKTYPE_STITCHED",
    "SDF_BLOCKTYPE_CONTIGUOUS",
    "SDF_BLOCKTYPE_LAGRANGIAN_MESH",
    "SDF_BLOCKTYPE_STATION",
    "SDF_BLOCKTYPE_STATION_DERIVED",
    "SDF_BLOCKTYPE_DATABLOCK",
    "SDF_BLOCKTYPE_NAMEVALUE",
};

const char *sdf_geometry_c[] = {
    "SDF_GEOMETRY_NULL",
    "SDF_GEOMETRY_CARTESIAN",
    "SDF_GEOMETRY_CYLINDRICAL",
    "SDF_GEOMETRY_SPHERICAL",
};

const char *sdf_stagger_c[] = {
    "SDF_STAGGER_CELL_CENTRE",
    "SDF_STAGGER_FACE_X",
    "SDF_STAGGER_FACE_Y",
    "SDF_STAGGER_EDGE_Z",
    "SDF_STAGGER_FACE_Z",
    "SDF_STAGGER_EDGE_Y",
    "SDF_STAGGER_EDGE_X",
    "SDF_STAGGER_VERTEX",
};

const char *sdf_datatype_c[] = {
    "SDF_DATATYPE_NULL",
    "SDF_DATATYPE_INTEGER4",
    "SDF_DATATYPE_INTEGER8",
    "SDF_DATATYPE_REAL4",
    "SDF_DATATYPE_REAL8",
    "SDF_DATATYPE_REAL16",
    "SDF_DATATYPE_CHARACTER",
    "SDF_DATATYPE_LOGICAL",
    "SDF_DATATYPE_OTHER",
};

const char *sdf_error_codes_c[] = {
    "SDF_ERR_SUCCESS",
    "SDF_ERR_ACCESS",
    "SDF_ERR_AMODE",
    "SDF_ERR_BAD_FILE",
    "SDF_ERR_CONVERSION",
    "SDF_ERR_DUP_DATAREP",
    "SDF_ERR_FILE",
    "SDF_ERR_FILE_EXISTS",
    "SDF_ERR_FILE_IN_USE",
    "SDF_ERR_INFO",
    "SDF_ERR_INFO_KEY",
    "SDF_ERR_INFO_NOKEY",
    "SDF_ERR_INFO_VALUE",
    "SDF_ERR_IO",
    "SDF_ERR_NOT_SAME",
    "SDF_ERR_NO_SPACE",
    "SDF_ERR_NO_SUCH_FILE",
    "SDF_ERR_QUOTA",
    "SDF_ERR_READ_ONLY",
    "SDF_ERR_UNSUPPORTED_DATAREP",
    "SDF_ERR_UNSUPPORTED_OPERATION",
    "SDF_ERR_UNKNOWN",
};

const int sdf_blocktype_len =
        sizeof(sdf_blocktype_c) / sizeof(sdf_blocktype_c[0]);
const int sdf_geometry_len =
        sizeof(sdf_geometry_c) / sizeof(sdf_geometry_c[0]);
const int sdf_stagger_len =
        sizeof(sdf_stagger_c) / sizeof(sdf_stagger_c[0]);
const int sdf_datatype_len =
        sizeof(sdf_datatype_c) / sizeof(sdf_datatype_c[0]);
const int sdf_error_codes_len =
        sizeof(sdf_error_codes_c) / sizeof(sdf_error_codes_c[0]);


static int sdf_abort(sdf_file_t *h)
{
#ifdef PARALLEL
    MPI_Abort(h->comm, 1);
#endif
    _exit(1);
    return 0;
}



int sdf_seek_set(sdf_file_t *h, off_t offset)
{
#ifdef PARALLEL
    return MPI_File_seek(h->filehandle, offset, MPI_SEEK_SET);
#else
    return FSEEKO(h->filehandle, offset, SEEK_SET);
#endif
}



int sdf_seek(sdf_file_t *h)
{
    return sdf_seek_set(h, h->current_location);
}



int sdf_broadcast(sdf_file_t *h, void *buf, int size)
{
#ifdef PARALLEL
    return MPI_Bcast(buf, size, MPI_BYTE, h->rank_master, h->comm);
#else
    return 0;
#endif
}



static int sdf_fopen(sdf_file_t *h, int mode)
{
    int ret = 0;

    // Abort for invalid mode argument
    assert(mode&SDF_READ || mode&SDF_WRITE);

#ifdef PARALLEL
    if (mode == SDF_READ)
        ret = MPI_File_open(h->comm, (char*)h->filename, MPI_MODE_RDONLY,
                MPI_INFO_NULL, &h->filehandle);
    else if (mode == SDF_WRITE)
        ret = MPI_File_open(h->comm, (char*)h->filename,
                MPI_MODE_WRONLY|MPI_MODE_CREATE, MPI_INFO_NULL, &h->filehandle);
    else if (mode == (SDF_READ|SDF_WRITE))
        ret = MPI_File_open(h->comm, (char*)h->filename,
                MPI_MODE_RDWR, MPI_INFO_NULL, &h->filehandle);
    else
        h->filehandle = 0;
#else
    if (mode == SDF_READ)
        h->filehandle = fopen(h->filename, "r");
    else if (mode == SDF_WRITE)
        h->filehandle = fopen(h->filename, "w");
    else if (mode == (SDF_READ|SDF_WRITE))
        h->filehandle = fopen(h->filename, "r+");
    else
        h->filehandle = NULL;
#endif
    if (!h->filehandle) ret = 1;

    return ret;
}



/** @ingroup control
 */
sdf_file_t *sdf_open(const char *filename, comm_t comm, int mode, int use_mmap)
{
    sdf_file_t *h;
    int ret;

    // Abort for invalid mode argument
    assert(mode&SDF_READ || mode&SDF_WRITE);

    // Create filehandle
    h = malloc(sizeof(*h));
    memset(h, 0, sizeof(*h));

#ifdef SDF_DEBUG
    h->dbg_count = DBG_CHUNK;
    h->dbg = h->dbg_buf = malloc(h->dbg_count);
#endif
    h->string_length = 64;
    h->indent = 0;

    h->done_header = 0;
    h->use_summary = 1;
    h->sdf_lib_version  = SDF_LIB_VERSION;
    h->sdf_lib_revision = SDF_LIB_REVISION;

    h->id_length = SDF_ID_LENGTH;
    // header length - must be updated if sdf_write_header changes
    h->first_block_location = SDF_HEADER_LENGTH;
    // block header length - must be updated if sdf_write_block_header changes
    h->block_header_length = SDF_BLOCK_HEADER_LENGTH;

#ifdef PARALLEL
    h->comm = comm;
    MPI_Comm_rank(h->comm, &h->rank);
    MPI_Comm_size(h->comm, &h->ncpus);
#else
    h->rank = 0;
    h->ncpus = 1;
#endif
    h->filename = malloc(strlen(filename)+1);
    memcpy(h->filename, filename, strlen(filename)+1);

    h->hashed_blocks_by_id = h->hashed_blocks_by_name = NULL;

    sdf_fopen(h, mode);
    if (!h->filehandle) {
        free(h->filename);
        free(h);
        h = NULL;
        return h;
    }

#ifndef PARALLEL
    if (use_mmap)
        h->mmap = "";
    else
#endif
        h->mmap = NULL;

    h->array_count = 20;

#ifndef PARALLEL
    if (h->mmap) {
        h->fd = fileno(h->filehandle);
    }
#endif

    if ((mode&SDF_READ) == 0) return h;

    ret = sdf_read_header(h);
    if (ret) {
        h = NULL;
        return h;
    }

    return h;
}



/** @ingroup control
 */
int sdf_close(sdf_file_t *h)
{
    sdf_extension_free_data(h);
    sdf_extension_unload();

    // No open file
    if (!h || !h->filehandle) return 1;

    sdf_fclose(h);

    // Destroy filehandle
    sdf_free_handle(h);

    return 0;
}



#define FREE_ITEM(value) do { \
        if ((value)) { \
            free((value)); \
            (value) = NULL; \
        } \
    } while(0)

#define FREE_ARRAY(block, value) do { \
    if (block->value) { \
        int _i; \
        for (_i = 0; _i < b->n##value; _i++) \
            if (block->value[_i]) free(block->value[_i]); \
        free(block->value); \
    }} while(0)



int sdf_free_block_data(sdf_file_t *h, sdf_block_t *b)
{
    int i;
    struct run_info *run;
    char **var;

    if (!b) return 1;

    if (b->grids) {
        if (!h->mmap && b->done_data && !b->dont_own_data) {
            for (i = 0; i < b->ngrids; i++) {
                if (b->grids[i]) {
                    if (b->data == b->grids[i])
                        b->data = NULL;
                    free(b->grids[i]);
                }
            }
        }
        free(b->grids);
        b->grids = NULL;
    }
    if (b->data && b->done_data && !b->dont_own_data) {
        if (b->blocktype == SDF_BLOCKTYPE_RUN_INFO) {
            run = b->data;
            if (run->commit_id)       free(run->commit_id);
            if (run->sha1sum)         free(run->sha1sum);
            if (run->compile_machine) free(run->compile_machine);
            if (run->compile_flags)   free(run->compile_flags);
            free(b->data);
            b->data = NULL;
        } else if (b->blocktype == SDF_BLOCKTYPE_NAMEVALUE) {
            if (b->datatype == SDF_DATATYPE_CHARACTER) {
                var = b->data;
                for (i=0; i < b->ndims; i++)
                    free(var[i]);
            }
            free(b->data);
            b->data = NULL;
        }
        if (!h->mmap && b->data) {
            free(b->data);
            b->data = NULL;
        }
    }
#if !defined(PARALLEL) && !defined(_WIN32)
    if (b->mmap) {
        munmap(b->mmap, b->mmap_len);
        b->mmap = NULL;
        b->mmap_len = 0;
    }
#endif
    if (b->node_list) free(b->node_list);
    if (b->boundary_cells) free(b->boundary_cells);
    b->node_list = NULL;
    b->boundary_cells = NULL;
    b->done_data = 0;
    b->grids = NULL;
    b->data = NULL;

    return 0;
}



int sdf_free_block(sdf_file_t *h, sdf_block_t *b)
{
    if (!b) return 1;

    sdf_delete_hash_block(h, b);

    FREE_ITEM(b->id);
    FREE_ITEM(b->units);
    FREE_ITEM(b->mesh_id);
    FREE_ITEM(b->material_id);
    FREE_ITEM(b->name);
    FREE_ITEM(b->material_name);
    FREE_ITEM(b->dim_mults);
    FREE_ITEM(b->extents);
    FREE_ITEM(b->station_nvars);
    FREE_ITEM(b->station_move);
    FREE_ITEM(b->station_x);
    FREE_ITEM(b->station_y);
    FREE_ITEM(b->station_z);
    FREE_ITEM(b->variable_types);
    FREE_ITEM(b->array_starts);
    FREE_ITEM(b->array_ends);
    FREE_ITEM(b->array_strides);
    FREE_ITEM(b->mimetype);
    FREE_ITEM(b->checksum_type);
    FREE_ITEM(b->checksum);
    FREE_ITEM(b->must_read);
    FREE_ITEM(b->station_id);
    FREE_ITEM(b->station_index);
    FREE_ITEM(b->obstacle_id);
    FREE_ITEM(b->vfm_id);

    FREE_ARRAY(b, station_ids);
    FREE_ARRAY(b, station_names);
    FREE_ARRAY(b, variable_ids);
    FREE_ARRAY(b, material_names);
    FREE_ARRAY(b, dim_labels);
    FREE_ARRAY(b, dim_units);

    sdf_free_block_data(h, b);

    free(b);
    b = NULL;

    return 0;
}



/** @ingroup control
 */
int sdf_free_blocklist_data(sdf_file_t *h)
{
    sdf_block_t *b, *next;
    int i;

    if (!h || !h->filehandle) return 1;

    // Destroy blocklist
    if (h->blocklist) {
        b = h->blocklist;
        for (i=0; i < h->nblocks; i++) {
            next = b->next;
            sdf_free_block_data(h, b);
            b = next;
        }
    }

    return 0;
}



static int sdf_free_handle(sdf_file_t *h)
{
    sdf_block_t *b, *next;
    int i;

    if (!h) return 1;

    // Destroy blocklist
    if (h->blocklist) {
        b = h->blocklist;
        for (i=0; i < h->nblocks; i++) {
            next = b->next;
            sdf_free_block(h, b);
            if (!next) break;
            b = next;
        }
    }
    // Destroy extension data
    if (h->ext_data) sdf_extension_free_data(h);
    // Destroy handle
    if (h->buffer) free(h->buffer);
    if (h->code_name) free(h->code_name);
    if (h->filename) free(h->filename);
    if (h->dbg_buf) free(h->dbg_buf);
    memset(h, 0, sizeof(sdf_file_t));
    free(h);
    h = NULL;

    return 0;
}



static int sdf_fclose(sdf_file_t *h)
{
    // No open file
    if (!h || !h->filehandle) return 1;

#ifdef PARALLEL
    MPI_Barrier(h->comm);

    MPI_File_close(&h->filehandle);
#else
    fclose(h->filehandle);
#endif
    h->filehandle = 0;

    return 0;
}



static int sdf_set_rank_master(sdf_file_t *h, int rank)
{
    if (h)
        h->rank_master = rank;
    else
        return -1;

    return 0;
}



static int sdf_read_nblocks(sdf_file_t *h)
{
    if (h)
        return h->nblocks;
    else
        return -1;
}



/*
static int sdf_read_jobid(sdf_file_t *h, sdf_jobid_t *jobid)
{
    if (h && jobid)
        memcpy(jobid, &h->jobid, sizeof(jobid));
    else
        return -1;

    return 0;
}
*/



#ifdef PARALLEL
static int factor2d(int ncpus, int64_t *dims, int *cpu_split)
{
    const int ndims = 2;
    int dmin[ndims], npoint_min[ndims], cpu_split_tmp[ndims], grids[ndims][2];
    int i, j, ii, jj, n, cpus, maxcpus, grid, split_big, dim;
    float gridav, deviation, mindeviation;

    cpus = 1;
    gridav = 1;
    for (i=0; i < ndims; i++) {
        dim = (int)dims[i];
        dmin[i] = MIN(ncpus, dim);
        cpus = cpus * dmin[i];
        gridav = gridav * dim;
    }
    mindeviation = gridav;
    gridav = gridav / ncpus;

    maxcpus = MIN(ncpus,cpus);

    for (j=0; j < dmin[1]; j++) {
        cpu_split_tmp[1] = dmin[1]-j;
    for (i=0; i < dmin[0]; i++) {
        cpu_split_tmp[0] = dmin[0]-i;

        cpus = 1;
        for (n=0; n < ndims; n++)
            cpus = cpus * cpu_split_tmp[n];

        if (cpus != maxcpus) continue;

        for (n=0; n < ndims; n++) {
            npoint_min[n] = (int)dims[n] / cpu_split_tmp[n];
            split_big = (int)dims[n] - cpu_split_tmp[n] * npoint_min[n];
            grids[n][0] = npoint_min[n];
            grids[n][1] = npoint_min[n] + 1;
            if (cpu_split_tmp[n] == split_big) grids[n][0] = 0;
            if (split_big == 0) grids[n][1] = 0;
        }

        for (ii=0; ii < 2; ii++) {
        for (jj=0; jj < 2; jj++) {
            grid = grids[0][ii] * grids[1][jj];
            deviation = ABS(grid-gridav);
            if (deviation < mindeviation) {
                mindeviation = deviation;
                for (n=0; n < ndims; n++)
                    cpu_split[n] = cpu_split_tmp[n];
            }
        }}
    }}

    return 0;
}



static int factor3d(int ncpus, int64_t *dims, int *cpu_split)
{
    const int ndims = 3;
    int dmin[ndims], npoint_min[ndims], cpu_split_tmp[ndims], grids[ndims][2];
    int i, j, k, ii, jj, kk, n, cpus, maxcpus, grid, split_big, dim;
    float gridav, deviation, mindeviation;

    cpus = 1;
    gridav = 1;
    for (i=0; i < ndims; i++) {
        dim = (int)dims[i];
        dmin[i] = MIN(ncpus, dim);
        cpus = cpus * dmin[i];
        gridav = gridav * dim;
    }
    mindeviation = gridav;
    gridav = gridav / ncpus;

    maxcpus = MIN(ncpus,cpus);

    for (k=0; k < dmin[2]; k++) {
        cpu_split_tmp[2] = dmin[2]-k;
    for (j=0; j < dmin[1]; j++) {
        cpu_split_tmp[1] = dmin[1]-j;
    for (i=0; i < dmin[0]; i++) {
        cpu_split_tmp[0] = dmin[0]-i;

        cpus = 1;
        for (n=0; n < ndims; n++)
            cpus = cpus * cpu_split_tmp[n];

        if (cpus != maxcpus) continue;

        for (n=0; n < ndims; n++) {
            npoint_min[n] = (int)dims[n] / cpu_split_tmp[n];
            split_big = (int)dims[n] - cpu_split_tmp[n] * npoint_min[n];
            grids[n][0] = npoint_min[n];
            grids[n][1] = npoint_min[n] + 1;
            if (cpu_split_tmp[n] == split_big) grids[n][0] = 0;
            if (split_big == 0) grids[n][1] = 0;
        }

        for (ii=0; ii < 2; ii++) {
        for (jj=0; jj < 2; jj++) {
        for (kk=0; kk < 2; kk++) {
            grid = grids[0][ii] * grids[1][jj] * grids[2][kk];
            deviation = ABS(grid-gridav);
            if (deviation < mindeviation) {
                mindeviation = deviation;
                for (n=0; n < ndims; n++)
                    cpu_split[n] = cpu_split_tmp[n];
            }
        }}}
    }}}

    return 0;
}
#endif



/** @ingroup control
 */
int sdf_get_domain_bounds(sdf_file_t *h, int rank, int *starts, int *local_dims)
{
    sdf_block_t *b = h->current_block;
    int n;
#ifdef PARALLEL
    int npoint_min, split_big, coords, div;
    int64_t old_dims[6];

    // Adjust dimensions to those of a cell-centred variable
    for (n = 0; n < b->ndims; n++) {
        old_dims[n] = b->dims[n];
        if (b->const_value[n]) b->dims[n]--;
        if (b->dims[n] < 1) b->dims[n] = 1;
    }

    memset(starts, 0, 3*sizeof(*starts));

    div = 1;
    for (n = 0; n < b->ndims; n++) {
        coords = (rank / div) % b->cpu_split[n];

        if (coords == 0)
            b->proc_min[n] = MPI_PROC_NULL;
        else
            b->proc_min[n] = rank - div;

        if (coords == b->cpu_split[n] - 1)
            b->proc_max[n] = MPI_PROC_NULL;
        else
            b->proc_max[n] = rank + div;

        div = div * b->cpu_split[n];
        npoint_min = (int)b->dims[n] / b->cpu_split[n];
        split_big = (int)(b->dims[n] - b->cpu_split[n] * npoint_min);
        if (coords >= split_big) {
            starts[n] = split_big * (npoint_min + 1)
                    + (coords - split_big) * npoint_min;
            local_dims[n] = npoint_min;
        } else {
            starts[n] = coords * (npoint_min + 1);
            local_dims[n] = npoint_min + 1;
        }

        // Add a layer of ghost cells for the VisIt reader
        if (h->internal_ghost_cells && !b->no_internal_ghost) {
            if (b->proc_min[n] != MPI_PROC_NULL) {
                b->ngb[2*n] = 1;
                local_dims[n]++;
                starts[n]--;
            }
            if (b->proc_max[n] != MPI_PROC_NULL) {
                b->ngb[2*n+1] = 1;
                local_dims[n]++;
            }
        }
    }

    // Return dimensions back to their original values
    for (n = 0; n < b->ndims; n++) {
        b->dims[n] = old_dims[n];
        // Add extra staggered value if required
        if (b->const_value[n]) local_dims[n]++;
    }
#else
    memset(starts, 0, 3*sizeof(*starts));
    for (n=0; n < b->ndims; n++) local_dims[n] = b->dims[n];
#endif
    for (n=b->ndims; n < 3; n++) local_dims[n] = 1;

    return 0;
}



int sdf_factor(sdf_file_t *h)
{
    sdf_block_t *b = h->current_block;
    int n;
#ifdef PARALLEL
    int64_t old_dims[6];
    int local_dims[SDF_MAXDIMS];

    // Adjust dimensions to those of a cell-centred variable
    for (n = 0; n < b->ndims; n++) {
        old_dims[n] = b->dims[n];
        if (b->const_value[n]) b->dims[n]--;
        if (b->dims[n] < 1) b->dims[n] = 1;
    }

    if (b->ndims == 2)
        factor2d(h->ncpus, b->dims, b->cpu_split);
    else
        factor3d(h->ncpus, b->dims, b->cpu_split);

    // Return dimensions back to their original values
    for (n = 0; n < b->ndims; n++)
        b->dims[n] = old_dims[n];

    sdf_get_domain_bounds(h, h->rank, b->starts, local_dims);
    for (n = 0; n < b->ndims; n++) {
        b->local_dims[n] = local_dims[n];
        if (b->local_dims[n] > b->dims[n])
            b->local_dims[n] = b->dims[n];
    }
    for (n = b->ndims; n < 3; n++)
        b->local_dims[n] = 1;
#else
    for (n = 0; n < 3; n++) b->local_dims[n] = b->dims[n];
#endif

    if (b->blocktype == SDF_BLOCKTYPE_PLAIN_MESH) {
        b->nelements_local = 0;
        for (n = 0; n < b->ndims; n++) b->nelements_local += b->local_dims[n];
    } else {
        b->nelements_local = 1;
        for (n = 0; n < b->ndims; n++) b->nelements_local *= b->local_dims[n];
    }

    return 0;
}



int sdf_convert_array_to_float(sdf_file_t *h, void **var_in, int count)
{
    sdf_block_t *b = h->current_block;

    if (h->swap) {
        if (b->datatype == SDF_DATATYPE_INTEGER4
                || b->datatype == SDF_DATATYPE_REAL4) {
            int i;
            int32_t *v = *var_in;
            for (i=0; i < count; i++) {
                _SDF_BYTE_SWAP32(*v);
                v++;
            }
        } else if (b->datatype == SDF_DATATYPE_INTEGER8
                || b->datatype == SDF_DATATYPE_REAL8) {
            int i;
            int64_t *v = *var_in;
            for (i=0; i < count; i++) {
                _SDF_BYTE_SWAP64(*v);
                v++;
            }
        }
    }

    if (h->use_float && b->datatype == SDF_DATATYPE_REAL8) {
        int i;
        float *r4;
        double *old_var, *r8;
        r8 = old_var = *var_in;
        r4 = *var_in = malloc(count * sizeof(float));
        for (i=0; i < count; i++)
            *r4++ = (float)(*r8++);
        if (!h->mmap) free(old_var);
        b->datatype_out = SDF_DATATYPE_REAL4;
#ifdef PARALLEL
        b->mpitype_out = MPI_FLOAT;
#endif
    }
    return 0;
}



int sdf_randomize_array(sdf_file_t *h, void **var_in, int count)
{
    sdf_block_t *b = h->current_block;

    sdf_random_init();

    if (b->datatype_out == SDF_DATATYPE_REAL8) {
        double tmp, *array = (double*)(*var_in);
        int i, id1, id2;

        for (i=0; i < count; i++) {
            id1 = 1LL * count * sdf_random() / SDF_MAX_RND;
            id2 = 1LL * count * sdf_random() / SDF_MAX_RND;
            tmp = array[id1];
            array[id1] = array[id2];
            array[id2] = tmp;
        }
    } else {
        float tmp, *array = (float*)(*var_in);
        int i, id1, id2;

        for (i=0; i < count; i++) {
            id1 = 1LL * count * sdf_random() / SDF_MAX_RND;
            id2 = 1LL * count * sdf_random() / SDF_MAX_RND;
            tmp = array[id1];
            array[id1] = array[id2];
            array[id2] = tmp;
        }
    }

    return 0;
}



static int sdf_header_copy(const sdf_file_t *h_in, sdf_file_t *h_out)
{
    sdf_file_t *h_tmp;

    if (h_in == h_out) return 1;

    h_tmp = malloc(sizeof(*h_tmp));
    memcpy(h_tmp, h_out, sizeof(*h_out));
    memcpy(h_out, h_in, sizeof(*h_in));

    h_out->dbg_count        = h_tmp->dbg_count;
    h_out->dbg              = h_tmp->dbg;
    h_out->dbg_buf          = h_tmp->dbg_buf;
    h_out->string_length    = h_tmp->string_length;
    h_out->indent           = h_tmp->indent;
    h_out->done_header      = h_tmp->done_header;
    h_out->use_summary      = h_tmp->use_summary;
    h_out->sdf_lib_version  = h_tmp->sdf_lib_version;
    h_out->sdf_lib_revision = h_tmp->sdf_lib_revision;
    h_out->comm             = h_tmp->comm;
    h_out->rank             = h_tmp->rank;
    h_out->ncpus            = h_tmp->ncpus;
    h_out->filename         = h_tmp->filename;
    h_out->filehandle       = h_tmp->filehandle;
    h_out->mmap             = h_tmp->mmap;

    free(h_tmp);

    return 0;
}



static uint32_t Q[41790], indx, carry, xcng, xs;

#define CNG (xcng = 69609 * xcng + 123)
#define XS (xs ^= xs<<13, xs ^= (unsigned)xs>>17, xs ^= xs>>5)
#define SUPR (indx < 41790 ? Q[indx++] : refill())
#define KISS SUPR + CNG + XS

static uint32_t refill(void)
{
    int i;
    uint64_t t;
    for (i=0; i < 41790; i++) {
        t = 7010176ULL * Q[i] + carry;
        carry = (t>>32);
        Q[i] = (uint32_t)~(t);
    }
    indx = 1;
    return (Q[0]);
}



static uint32_t sdf_random(void)
{
    return KISS;
}



static void sdf_random_init(void)
{
    int i;
    indx = 41790;
    carry = 362436;
    xcng = 1236789;
    xs = 521288629;
    for (i=0; i < 41790; i++) Q[i] = CNG + XS;
    for (i=0; i < 41790; i++) sdf_random();
}



// Currently only works in serial
int sdf_block_set_array_section(sdf_block_t *b, const int ndims,
                                const int64_t *starts, const int64_t *ends,
                                const int64_t *strides)
{
    int64_t nelements_local;
    int i, ndims_min;

    if (b->ndims < 1) return 1;

    if (b->blocktype == SDF_BLOCKTYPE_PLAIN_MESH
            || b->blocktype == SDF_BLOCKTYPE_CPU_SPLIT) {
        nelements_local = 0;
        for (i = 0; i < b->ndims; i++)
            nelements_local += b->local_dims[i];
    } else if (b->blocktype == SDF_BLOCKTYPE_POINT_MESH) {
        nelements_local = b->local_dims[0];
    } else {
        nelements_local = 1;
        for (i = 0; i < b->ndims; i++)
            nelements_local *= b->local_dims[i];
    }

    if (!starts && !ends) {
        b->nelements_local = nelements_local;
        FREE_ITEM(b->array_starts);
        FREE_ITEM(b->array_ends);
        FREE_ITEM(b->array_strides);

        return 0;
    }

    if (!b->array_starts) {
        b->array_starts  = calloc(b->ndims, sizeof(*b->array_starts));
        b->array_ends    = malloc(b->ndims * sizeof(*b->array_ends));
        b->array_strides = malloc(b->ndims * sizeof(*b->array_strides));

        for (i = 0; i < b->ndims; i++) {
            b->array_ends[i] = b->local_dims[i];
            b->array_strides[i] = 1;
        }
    }

    ndims_min = (ndims < b->ndims) ? ndims : b->ndims;

    if (starts) {
        for (i = 0; i < ndims_min; i++) {
            if (starts[i] < 0)
                b->array_starts[i] = starts[i] + b->local_dims[i];
            else if (starts[i] > b->local_dims[i])
                b->array_starts[i] = b->local_dims[i];
            else
                b->array_starts[i] = starts[i];
        }
    }

    if (ends) {
        for (i = 0; i < ndims_min; i++) {
            if (ends[i] < 0)
                b->array_ends[i] = ends[i] + b->local_dims[i];
            else if (ends[i] == 0 && starts && starts[i] == -1)
                b->array_ends[i] = b->local_dims[i];
            else if (ends[i] < b->local_dims[i])
                b->array_ends[i] = ends[i];
        }
    } else {
        for (i = 0; i < ndims_min; i++)
            b->array_ends[i] = b->array_starts[i] + 1;
    }

    if (strides) {
        for (i = 0; i < ndims_min; i++) {
            if (strides[i] != 1) {
                fprintf(stderr,
                        "Array section striding is not yet implemented\n");
                break;
            }
        }
    }

    if (b->blocktype == SDF_BLOCKTYPE_PLAIN_MESH
            || b->blocktype == SDF_BLOCKTYPE_POINT_MESH) {
        b->nelements_local = 0;
        for (i = 0; i < b->ndims; i++)
            b->nelements_local += (b->array_ends[i] - b->array_starts[i]);
    } else {
        b->nelements_local = 1;
        for (i = 0; i < b->ndims; i++)
            b->nelements_local *= (b->array_ends[i] - b->array_starts[i]);
    }

    if (b->nelements_local == nelements_local) {
        FREE_ITEM(b->array_starts);
        FREE_ITEM(b->array_ends);
        FREE_ITEM(b->array_strides);
    }

    return 0;
}


char *sdf_get_library_commit_id(void)
{
    return SDF_COMMIT_ID;
}


char *sdf_get_library_commit_date(void)
{
    return SDF_COMMIT_DATE;
}


int sdf_has_debug_info(void)
{
#ifdef SDF_DEBUG
    return 1;
#else
    return 0;
#endif
}
