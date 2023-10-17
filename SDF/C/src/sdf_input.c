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
#include "sdf_input_point.h"
#include "sdf_input_station.h"
#include "sdf_control.h"

#ifdef PARALLEL
# include <mpi.h>
#else
# ifdef _WIN32
#  include <io.h>
# else
#  include <unistd.h>
#  include <sys/mman.h>
# endif
#endif

/**
 * @defgroup input
 * @brief Routines for reading from an SDF file
 */

#define SDF_COMMON_INFO() do { \
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


static int sdf_read_constant(sdf_file_t *h);
static int sdf_read_stitched(sdf_file_t *h);
static int sdf_read_stitched_material(sdf_file_t *h);
static int sdf_read_stitched_matvar(sdf_file_t *h);
static int sdf_read_stitched_species(sdf_file_t *h);
static int sdf_read_stitched_obstacle_group(sdf_file_t *h);
static int sdf_read_array(sdf_file_t *h);
static int sdf_read_array_info(sdf_file_t *h);
static int sdf_read_cpu_split_info(sdf_file_t *h);
static int sdf_read_run_info(sdf_file_t *h);
static int sdf_read_datablock_info(sdf_file_t *h);
static int sdf_read_datablock(sdf_file_t *h);
static int sdf_read_namevalue(sdf_file_t *h);
static void build_summary_buffer(sdf_file_t *h);
static int sdf_read_next_block_header(sdf_file_t *h);



void sdf_trim(char *str)
{
    int i, len = strlen(str);
    char *ptr = str + len - 1;

    for (i=0, ptr=str+len-1; i < len && *ptr == ' '; i++, ptr--)
        *ptr = '\0';
}



static inline int sdf_get_next_block(sdf_file_t *h)
{
    if (h->blocklist) {
        if (!h->current_block)
            h->current_block = h->blocklist;
        else if (h->current_block->next)
            h->current_block = h->current_block->next;
        else {
            sdf_block_t *block = malloc(sizeof(sdf_block_t));
            memset(block, 0, sizeof(sdf_block_t));
            block->block_start = h->tail->next_block_location;
            if (h->use_summary)
                block->summary_block_start = block->block_start;
            else
                block->inline_block_start = block->block_start;
            h->tail->next = block;
            h->tail->next->prev = h->tail;
            h->current_block = h->tail = block;
        }
    } else {
        sdf_block_t *block = malloc(sizeof(sdf_block_t));
        memset(block, 0, sizeof(sdf_block_t));
        if (h->use_summary)
            block->summary_block_start = block->block_start
                    = h->summary_location;
        else
            block->inline_block_start = block->block_start
                    = h->first_block_location;
        h->blocklist = h->tail = h->current_block = block;
    }

    h->current_block->in_file = 1;

    return 0;
}



int sdf_read_bytes(sdf_file_t *h, char *buf, size_t buflen)
{
#ifdef PARALLEL
    return MPI_File_read(h->filehandle, buf, (int)buflen, MPI_BYTE,
            MPI_STATUS_IGNORE);
#else
    return (1 != fread(buf, buflen, 1, h->filehandle));
#endif
}



/** @ingroup input
 *  @{
 */
int sdf_read_header(sdf_file_t *h)
{
    size_t buflen;

    h->indent = 0;

    if (h->done_header) return 1;

    buflen = h->first_block_location;
    h->buffer = malloc(buflen);

    h->current_location = h->start_location = 0;

    if (h->rank == h->rank_master) {
        sdf_seek(h);
        sdf_read_bytes(h, h->buffer, buflen);
    }
    sdf_broadcast(h, h->buffer, buflen);

    // If this isn't SDF_MAGIC then this isn't a SDF file;
    if (memcmp(h->buffer, SDF_MAGIC, 4) != 0) {
        sdf_close(h);
        return -1;
    }

    h->current_location += 4;

    SDF_READ_ENTRY_INT4(h->endianness);
    if (h->endianness == 0x0f0e0201) h->swap = 1;

    SDF_READ_ENTRY_INT4(h->file_version);
    if (h->file_version > SDF_VERSION) {
        sdf_close(h);
        return -1;
    }

    SDF_READ_ENTRY_INT4(h->file_revision);

    SDF_READ_ENTRY_ID(h->code_name);

    SDF_READ_ENTRY_INT8(h->first_block_location);

    SDF_READ_ENTRY_INT8(h->summary_location);

    SDF_READ_ENTRY_INT4(h->summary_size);

    SDF_READ_ENTRY_INT4(h->nblocks_file);

    SDF_READ_ENTRY_INT4(h->block_header_length);

    SDF_READ_ENTRY_INT4(h->step);

    SDF_READ_ENTRY_REAL8(h->time);

    SDF_READ_ENTRY_INT4(h->jobid1);

    SDF_READ_ENTRY_INT4(h->jobid2);

    SDF_READ_ENTRY_INT4(h->string_length);

    SDF_READ_ENTRY_INT4(h->code_io_version);

    SDF_READ_ENTRY_LOGICAL(h->restart_flag);

    SDF_READ_ENTRY_LOGICAL(h->other_domains);

    SDF_READ_ENTRY_LOGICAL(h->station_file);

    free(h->buffer);
    h->buffer = NULL;

    h->current_location = h->first_block_location;
    h->done_header = 1;
    h->nblocks = h->nblocks_file;

    if (h->summary_size == 0) {
        h->use_summary = 0;
        h->summary_metadata_invalid = 1;
    }

    return 0;
}



int sdf_read_summary(sdf_file_t *h)
{
    if (h->blocklist && !h->tmp_flag) {
        h->current_block = h->blocklist;
        return 0;
    }

    if (!h->done_header) sdf_read_header(h);
    h->current_block = NULL;

    // Read the whole summary block into a temporary buffer on rank 0
    if (h->use_summary > 0) {
        h->current_location = h->start_location = h->summary_location;
        h->buffer = malloc(h->summary_size);

        if (h->rank == h->rank_master) {
            sdf_seek(h);
            sdf_read_bytes(h, h->buffer, h->summary_size);
        }
        h->summary_metadata_read = 1;
    } else {
        build_summary_buffer(h);
        h->current_location = h->start_location = h->first_block_location;
    }

    // Send the temporary buffer to all processors
    sdf_broadcast(h, h->buffer, h->summary_size);

    return 0;
}



int sdf_purge_duplicates(sdf_file_t *h)
{
    sdf_block_t *b, *next, *subnext;
    char *id;
    int pos, i;

    // Ensure all blocks have been hashed first
    sdf_hash_block_list(h);

    // Remove or mangle the IDs of any duplicated blocks
    next = h->blocklist;
    while (next) {
        b = next;
        next = b->next;

        subnext = sdf_find_block_by_id(h, b->id);
        if (!subnext || subnext == b)
            continue;

        sdf_delete_hash_block(h, b);

        if (h->purge_duplicated_ids) {
            sdf_modify_remove_block(h, b);
        } else {
            pos = strlen(b->id);
            if (pos == h->id_length)
                pos--;

            id = calloc(pos + 3, 1);
            strncpy(id, b->id, pos+3);
            free(b->id);
            b->id = id;

            // Mangle ID by appending an integer. Allow for up to 99
            // duplicated blocks
            for (i=1; i < 99; i++) {
                if (i == 10 && pos == h->id_length - 1)
                    pos--;

                sprintf(b->id + pos, "%d", i);

                // Check that the new ID is unique
                subnext = sdf_find_block_by_id(h, b->id);
                if (!subnext)
                    break;
            }

            if (subnext)
                // Discard block if we can't find a unique ID
                sdf_modify_remove_block(h, b);
            else
                sdf_hash_block(h, b);
        }
    }

    return 0;
}



int sdf_read_blocklist(sdf_file_t *h)
{
    int i;
#ifdef PARALLEL
    int fix;
    sdf_block_t *b, *next, *mesh;
#endif

    sdf_read_summary(h);

    // Construct the metadata blocklist using the contents of the buffer
    for (i = 0; i < h->nblocks; i++) {
        SDF_DPRNT("\n");
        sdf_read_block_info(h);
    }

    free(h->buffer);

    h->buffer = NULL;
    h->current_block = h->blocklist;
    h->last_block_in_file = h->tail;

    sdf_purge_duplicates(h);

#ifdef PARALLEL
    // Hack to fix cartesian blocks whose mesh sizes don't match the stagger
    next = h->blocklist;
    while (next) {
        b = next;
        next = b->next;
        if (b->blocktype == SDF_BLOCKTYPE_PLAIN_VARIABLE && b->stagger) {
            fix = 0;
            mesh = sdf_find_block_by_id(h, b->mesh_id);
            for (i = 0; i < b->ndims; i++) {
                if (b->const_value[i] && b->dims[i] != mesh->dims[i]) {
                    fix = 1;
                    b->const_value[i] = 0;
                }
            }
            if (fix) {
                // Re-calculate per block parallel factorisation
                sdf_factor(h);
            }
        }
    }
#endif
    return 0;
}



int sdf_read_block_info(sdf_file_t *h)
{
    sdf_block_t *b;
    int ret = 0;

    sdf_read_next_block_header(h);
    b = h->current_block;
    if (b->done_info) return 0;

    h->indent += 2;
    if (b->blocktype == SDF_BLOCKTYPE_PLAIN_MESH
            || b->blocktype == SDF_BLOCKTYPE_LAGRANGIAN_MESH)
        ret = sdf_read_plain_mesh_info(h);
    else if (b->blocktype == SDF_BLOCKTYPE_POINT_MESH)
        ret = sdf_read_point_mesh_info(h);
    else if (b->blocktype == SDF_BLOCKTYPE_PLAIN_VARIABLE)
        ret = sdf_read_plain_variable_info(h);
    else if (b->blocktype == SDF_BLOCKTYPE_POINT_VARIABLE)
        ret = sdf_read_point_variable_info(h);
    else if (b->blocktype == SDF_BLOCKTYPE_CONSTANT)
        ret = sdf_read_constant(h);
    else if (b->blocktype == SDF_BLOCKTYPE_ARRAY)
        ret = sdf_read_array_info(h);
    else if (b->blocktype == SDF_BLOCKTYPE_CPU_SPLIT)
        ret = sdf_read_cpu_split_info(h);
    else if (b->blocktype == SDF_BLOCKTYPE_RUN_INFO)
        ret = sdf_read_run_info(h);
    else if (b->blocktype == SDF_BLOCKTYPE_STITCHED
            || b->blocktype == SDF_BLOCKTYPE_CONTIGUOUS
            || b->blocktype == SDF_BLOCKTYPE_STITCHED_TENSOR
            || b->blocktype == SDF_BLOCKTYPE_CONTIGUOUS_TENSOR)
        ret = sdf_read_stitched(h);
    else if (b->blocktype == SDF_BLOCKTYPE_STITCHED_MATERIAL
            || b->blocktype == SDF_BLOCKTYPE_CONTIGUOUS_MATERIAL)
        ret = sdf_read_stitched_material(h);
    else if (b->blocktype == SDF_BLOCKTYPE_STITCHED_MATVAR
            || b->blocktype == SDF_BLOCKTYPE_CONTIGUOUS_MATVAR)
        ret = sdf_read_stitched_matvar(h);
    else if (b->blocktype == SDF_BLOCKTYPE_STITCHED_SPECIES
            || b->blocktype == SDF_BLOCKTYPE_CONTIGUOUS_SPECIES)
        ret = sdf_read_stitched_species(h);
    else if (b->blocktype == SDF_BLOCKTYPE_STITCHED_OBSTACLE_GROUP)
        ret = sdf_read_stitched_obstacle_group(h);
    else if (b->blocktype == SDF_BLOCKTYPE_STATION)
        ret = sdf_read_station_info(h);
    else if (b->blocktype == SDF_BLOCKTYPE_DATABLOCK)
        ret = sdf_read_datablock_info(h);
    else if (b->blocktype == SDF_BLOCKTYPE_NAMEVALUE)
        ret = sdf_read_namevalue(h);

    // Fix up block_start values for inline metadata
    if (!h->use_summary) {
        if (b->prev)
            b->block_start = b->prev->next_block_location;
        else if (b == h->blocklist)
            b->block_start = h->first_block_location;
    }

    return ret;
}



int sdf_read_data(sdf_file_t *h)
{
    sdf_block_t *b;

    b = h->current_block;

    if (b->populate_data) {
        b->populate_data(h, b);
        return 0;
    } else if (b->blocktype == SDF_BLOCKTYPE_PLAIN_MESH)
        return sdf_read_plain_mesh(h);
    else if (b->blocktype == SDF_BLOCKTYPE_LAGRANGIAN_MESH)
        return sdf_read_lagran_mesh(h);
    else if (b->blocktype == SDF_BLOCKTYPE_POINT_MESH)
        return sdf_read_point_mesh(h);
    else if (b->blocktype == SDF_BLOCKTYPE_PLAIN_VARIABLE)
        return sdf_read_plain_variable(h);
    else if (b->blocktype == SDF_BLOCKTYPE_POINT_VARIABLE)
        return sdf_read_point_variable(h);
    else if (b->blocktype == SDF_BLOCKTYPE_ARRAY
            || b->blocktype == SDF_BLOCKTYPE_CPU_SPLIT)
        return sdf_read_array(h);
    else if (b->blocktype == SDF_BLOCKTYPE_DATABLOCK)
        return sdf_read_datablock(h);

    return 1;
}
/** @} */



// Read the block header into the current block
static int sdf_read_next_block_header(sdf_file_t *h)
{
    sdf_block_t *b;
    int i;

    if (!h->done_header) {
        if (h->rank == h->rank_master) {
            fprintf(stderr, "*** ERROR ***\n");
            fprintf(stderr, "SDF header has not been read. Ignoring call.\n");
        }
        return 1;
    }

    sdf_get_next_block(h);
    b = h->current_block;

    if (!b) {
        if (h->rank == h->rank_master) {
            fprintf(stderr, "*** ERROR ***\n");
            fprintf(stderr, "SDF block not initialised. Ignoring call.\n");
        }
        return 1;
    }

    if (b->done_header) {
        h->current_location = b->block_start + h->block_header_length;
        return 0;
    }

    h->indent = 2;

    h->current_location = b->block_start;

    SDF_READ_ENTRY_INT8(b->next_block_location);

    SDF_READ_ENTRY_INT8(b->data_location);

    SDF_READ_ENTRY_ID(b->id);

    SDF_READ_ENTRY_INT8(b->data_length);

    SDF_READ_ENTRY_TYPE(blocktype);

    SDF_READ_ENTRY_TYPE(datatype);

    SDF_READ_ENTRY_INT4(b->ndims);

    SDF_READ_ENTRY_STRING(b->name);

    // Older versions of the file did not contain the block
    // info length in the header.
    if (h->file_version + h->file_revision > 1)
        SDF_READ_ENTRY_INT4(b->info_length);

    if (h->use_summary) {
        b->summary_next_block_location = b->next_block_location;
    } else {
        b->inline_next_block_location = b->next_block_location;
        b->next_block_location = b->block_start
                + h->block_header_length + b->info_length;
    }

    if (b->blocktype == SDF_BLOCKTYPE_POINT_VARIABLE
            || b->blocktype == SDF_BLOCKTYPE_POINT_MESH)
        b->stagger = SDF_STAGGER_VERTEX;
    else
        b->stagger = SDF_STAGGER_CELL_CENTRE;
    for (i = 0; i < 3; i++) b->dims[i] = 1;

    if (b->blocktype == SDF_BLOCKTYPE_STATION)
        h->station_file = 1;

    b->done_header = 1;
    h->current_location = b->block_start + h->block_header_length;

    b->datatype_out = b->datatype;
    if (h->use_float && b->datatype == SDF_DATATYPE_REAL8)
        b->datatype_out = SDF_DATATYPE_REAL4;
#ifdef PARALLEL
    switch (b->datatype) {
    case(SDF_DATATYPE_REAL4):
        b->mpitype = MPI_FLOAT;
        break;
    case(SDF_DATATYPE_REAL8):
        b->mpitype = MPI_DOUBLE;
        break;
    case(SDF_DATATYPE_INTEGER4):
        b->mpitype = MPI_INT;
        break;
    case(SDF_DATATYPE_INTEGER8):
        b->mpitype = MPI_LONG_LONG;
        break;
    case(SDF_DATATYPE_CHARACTER):
        b->mpitype = MPI_CHAR;
        break;
    case(SDF_DATATYPE_LOGICAL):
        b->mpitype = MPI_CHAR;
        break;
    }
    b->mpitype_out = b->mpitype;
#endif

    sdf_hash_block(h, b);

    return 0;
}



// Read all the block metadata sections and copy them into
// one contiguous buffer (h->buffer).
static void build_summary_buffer(sdf_file_t *h)
{
    int count, buflen;
    int64_t data_location, block_location, next_block_location;
    int32_t info_length, nblocks;
    char *bufptr;
    char skip_summary;

    struct list_entry {
        void *buffer;
        int len;
        struct list_entry *next;
    } *blockbuf_head, *blockbuf;

    h->current_location = h->first_block_location;

    // Read the file and build the buffer on rank zero.
    if (h->rank == h->rank_master) {
        next_block_location = block_location = h->current_location;

        blockbuf_head = blockbuf = calloc(1,sizeof(*blockbuf));
        buflen = 0;
        nblocks = 0;
        skip_summary = (h->summary_location
               && h->summary_location != h->first_block_location);
        // Read the block metadata into a temporary linked list structure
        while (1) {
            if (skip_summary
                    && h->current_location >= h->summary_location) break;

            sdf_seek(h);

            // Read the fixed length block header
            blockbuf->len = h->block_header_length;
            blockbuf->buffer = malloc(blockbuf->len);
            count = sdf_read_bytes(h, blockbuf->buffer, blockbuf->len);
            if (count != 0) break;

            memcpy(&next_block_location, blockbuf->buffer,
                    sizeof(next_block_location));
            memcpy(&data_location, (char*)blockbuf->buffer + SOI8,
                    sizeof(data_location));

            if (h->swap) {
                _SDF_BYTE_SWAP64(next_block_location);
                _SDF_BYTE_SWAP64(data_location);
            }

            // Older versions of the file did not contain the block
            // info length in the header.
            if (h->file_version + h->file_revision > 1) {
                memcpy(&info_length,
                        (char*)blockbuf->buffer + 3 * (SOI4 + SOI8)
                        + h->string_length + h->id_length,
                        sizeof(info_length));
                if (h->swap) _SDF_BYTE_SWAP32(info_length);
            } else {
                if (data_location > block_location)
                    info_length = (int32_t)(data_location
                            - block_location) - h->block_header_length;
                else
                    info_length = (int32_t)(next_block_location
                            - block_location) - h->block_header_length;
            }

            // Read the block specific metadata if it exists
            if (info_length > 0) {
                blockbuf->next = calloc(1,sizeof(*blockbuf));
                blockbuf = blockbuf->next;

                blockbuf->len = info_length;
                blockbuf->buffer = malloc(blockbuf->len);
                count = sdf_read_bytes(h, blockbuf->buffer, blockbuf->len);
                if (count != 0) break;
            }

            buflen += h->block_header_length + info_length;

            nblocks++;

            blockbuf->next = calloc(1,sizeof(*blockbuf));
            blockbuf = blockbuf->next;

            if (h->current_location > next_block_location) break;
            if (h->ignore_nblocks == 0 && nblocks >= h->nblocks) break;

            h->current_location = block_location = next_block_location;
        }

        if (h->ignore_nblocks != 0) h->nblocks = nblocks;

        if (blockbuf->buffer) {
            free(blockbuf->buffer);
            blockbuf->buffer = NULL;
        }

        // Copy the contents of the linked list into a single contiguous buffer.
        bufptr = h->buffer = malloc(buflen);
        while (blockbuf_head) {
            blockbuf = blockbuf_head;

            if (blockbuf->buffer == NULL) {
                free(blockbuf);
                break;
            }

            memcpy(bufptr, blockbuf->buffer, blockbuf->len);
            bufptr += blockbuf->len;

            blockbuf_head = blockbuf->next;

            free(blockbuf->buffer);
            free(blockbuf);
        }
    }

    // Send the temporary buffer length to all processors
    sdf_broadcast(h, &buflen, sizeof(buflen));

    // Allocate the buffer on all non rank zero processors
    if (h->rank != h->rank_master)
        h->buffer = malloc(buflen);

    h->summary_size = buflen;
    h->current_location = 0;
    h->inline_metadata_read = 1;
}



static int sdf_read_stitched(sdf_file_t *h)
{
    sdf_block_t *b;

    SDF_COMMON_INFO();

    // Metadata is
    // - stagger   INTEGER(i4)
    // - meshid    CHARACTER(id_length)
    // - varids    ndims*CHARACTER(id_length)

    SDF_READ_ENTRY_INT4(b->stagger);

    SDF_READ_ENTRY_ID(b->mesh_id);

    SDF_READ_ENTRY_ARRAY_ID(b->variable_ids, b->ndims);
    b->nvariable_ids = b->ndims;

    if (b->blocktype == SDF_BLOCKTYPE_STITCHED
            || b->blocktype == SDF_BLOCKTYPE_STITCHED_TENSOR) b->done_data = 1;

    return 0;
}



static int sdf_read_stitched_material(sdf_file_t *h)
{
    sdf_block_t *b;

    SDF_COMMON_INFO();

    // Metadata is
    // - stagger   INTEGER(i4)
    // - meshid    CHARACTER(id_length)
    // - matnames  ndims*CHARACTER(string_length)
    // - varids    ndims*CHARACTER(id_length)

    SDF_READ_ENTRY_INT4(b->stagger);

    SDF_READ_ENTRY_ID(b->mesh_id);

    SDF_READ_ENTRY_ARRAY_STRING(b->material_names, b->ndims);

    SDF_READ_ENTRY_ARRAY_ID(b->variable_ids, b->ndims);
    b->nmaterial_names = b->nvariable_ids = b->ndims;

    b->done_data = 1;

    return 0;
}



static int sdf_read_stitched_matvar(sdf_file_t *h)
{
    sdf_block_t *b;

    SDF_COMMON_INFO();

    // Metadata is
    // - stagger   INTEGER(i4)
    // - meshid    CHARACTER(id_length)
    // - matid     CHARACTER(id_length)
    // - varids    ndims*CHARACTER(id_length)

    SDF_READ_ENTRY_INT4(b->stagger);

    SDF_READ_ENTRY_ID(b->mesh_id);

    SDF_READ_ENTRY_ID(b->material_id);

    SDF_READ_ENTRY_ARRAY_ID(b->variable_ids, b->ndims);
    b->nvariable_ids = b->ndims;

    b->done_data = 1;

    return 0;
}



static int sdf_read_stitched_species(sdf_file_t *h)
{
    sdf_block_t *b;

    SDF_COMMON_INFO();

    // Metadata is
    // - stagger   INTEGER(i4)
    // - meshid    CHARACTER(id_length)
    // - matid     CHARACTER(id_length)
    // - matname   CHARACTER(string_length)
    // - specnames ndims*CHARACTER(string_length)
    // - varids    ndims*CHARACTER(id_length)

    SDF_READ_ENTRY_INT4(b->stagger);

    SDF_READ_ENTRY_ID(b->mesh_id);

    SDF_READ_ENTRY_ID(b->material_id);

    SDF_READ_ENTRY_STRING(b->material_name);

    SDF_READ_ENTRY_ARRAY_STRING(b->material_names, b->ndims);

    SDF_READ_ENTRY_ARRAY_ID(b->variable_ids, b->ndims);
    b->nvariable_ids = b->nmaterial_names = b->ndims;

    b->done_data = 1;

    return 0;
}



static int sdf_read_stitched_obstacle_group(sdf_file_t *h)
{
    sdf_block_t *b;

    SDF_COMMON_INFO();

    // Metadata is
    // - stagger         INTEGER(i4)
    // - obstacle_id     CHARACTER(id_length)
    // - vfm_id          CHARACTER(id_length)
    // - obstacle_names  ndims*CHARACTER(string_length)

    SDF_READ_ENTRY_INT4(b->stagger);

    SDF_READ_ENTRY_ID(b->obstacle_id);

    SDF_READ_ENTRY_ID(b->vfm_id);

    SDF_READ_ENTRY_ARRAY_STRING(b->material_names, b->ndims);
    b->nmaterial_names = b->ndims;

    b->done_data = 1;

    return 0;
}



static int sdf_read_constant(sdf_file_t *h)
{
    sdf_block_t *b;

    // Metadata is
    // - value     TYPE_SIZE

    b = h->current_block;
    SDF_READ_ENTRY_CONST(b->const_value);

    b->stagger = SDF_STAGGER_VERTEX;

    return 0;
}



static int sdf_array_datatype(sdf_file_t *h)
{
    sdf_block_t *b = h->current_block;
    int n;

    for (n=0; n < b->ndims; n++) b->local_dims[n] = b->dims[n];
    for (n=b->ndims; n < 3; n++) b->local_dims[n] = 1;

    b->nelements_local = 1;
    for (n=0; n < b->ndims; n++) b->nelements_local *= b->local_dims[n];
#ifdef PARALLEL
    MPI_Type_contiguous(b->nelements_local, b->mpitype, &b->distribution);
    MPI_Type_commit(&b->distribution);
#endif

    return 0;
}



static int sdf_read_run_info(sdf_file_t *h)
{
    sdf_block_t *b;
    struct run_info *run;

    run = calloc(1, sizeof(*run));

    // Metadata is
    // - version   INTEGER(i4)
    // - revision  INTEGER(i4)
    // - commit_id CHARACTER(string_length)
    // - sha1sum   CHARACTER(string_length)
    // - compmac   CHARACTER(string_length)
    // - compflag  CHARACTER(string_length)
    // - defines   INTEGER(i8)
    // - compdate  INTEGER(i4)
    // - rundate   INTEGER(i4)
    // - iodate    INTEGER(i4)
    // - minor_rev INTEGER(i4)

    SDF_COMMON_INFO();
    SDF_READ_ENTRY_INT4(run->version);
    SDF_READ_ENTRY_INT4(run->revision);
    SDF_READ_ENTRY_STRING(run->commit_id);
    SDF_READ_ENTRY_STRING(run->sha1sum);
    SDF_READ_ENTRY_STRING(run->compile_machine);
    SDF_READ_ENTRY_STRING(run->compile_flags);
    SDF_READ_ENTRY_INT8(run->defines);
    SDF_READ_ENTRY_INT4(run->compile_date);
    SDF_READ_ENTRY_INT4(run->run_date);
    SDF_READ_ENTRY_INT4(run->io_date);
    if (h->file_version == 1 && h->file_revision < 2)
        run->minor_rev = 0;
    else
        SDF_READ_ENTRY_INT4(run->minor_rev);

    h->current_block->data = run;
    b->done_data = 1;

    return 0;
}



static int sdf_read_array_info(sdf_file_t *h)
{
    sdf_block_t *b;
    int i;
    int32_t dims_in[SDF_MAXDIMS];
    int32_t *dims_ptr = dims_in;

    // Metadata is
    // - dims      INTEGER(i4), DIMENSION(ndims)

    SDF_COMMON_INFO();

    SDF_READ_ENTRY_ARRAY_INT4(dims_ptr, b->ndims);
    b->nelements = 1;
    for (i = 0; i < b->ndims; i++) {
        b->local_dims[i] = b->dims[i] = dims_in[i];
        b->nelements *= b->dims[i];
    }
    b->nelements_local = b->nelements;

    return 0;
}



static int sdf_read_cpu_split_info(sdf_file_t *h)
{
    sdf_block_t *b;
    int i;
    int32_t dims_in[SDF_MAXDIMS];
    int32_t *dims_ptr = dims_in;

    // Metadata is
    // - dims      INTEGER(i4), DIMENSION(ndims)

    SDF_COMMON_INFO();

    SDF_READ_ENTRY_INT4(b->geometry);
    SDF_READ_ENTRY_ARRAY_INT4(dims_ptr, b->ndims);
    for (i = 0; i < b->ndims; i++) {
        b->local_dims[i] = b->dims[i] = dims_in[i];
    }
    if (b->geometry == 1 || b->geometry == 4) {
        b->nelements_local = 0;
        for (i = 0; i < b->ndims; i++)
            b->nelements_local += b->dims[i];
    } else if (b->geometry == 2) {
        b->nelements_local = (int)(b->dims[0] * (b->dims[1] + 1));
        if (b->ndims > 2)
            b->nelements_local += b->dims[0] * b->dims[1] * b->dims[2];
    } else if (b->geometry == 3) {
        b->nelements_local = 1;
        for (i = 0; i < b->ndims; i++)
            b->nelements_local *= b->dims[i];
    }

    return 0;
}



static int sdf_read_datablock_info(sdf_file_t *h)
{
    sdf_block_t *b;

    // Metadata is
    // - mimetype       CHARACTER(id_length)
    // - checksum_type  CHARACTER(id_length)
    // - checksum       CHARACTER(string_length)

    SDF_COMMON_INFO();

    SDF_READ_ENTRY_ID(b->mimetype);
    SDF_READ_ENTRY_ID(b->checksum_type);
    SDF_READ_ENTRY_STRING(b->checksum);

    return 0;
}



static int sdf_read_namevalue(sdf_file_t *h)
{
    sdf_block_t *b;
    int32_t *i4 = NULL;
    int64_t *i8 = NULL;
    float *r4 = NULL;
    double *r8 = NULL;
    char *logical = NULL;
    char **string = NULL;

    // Metadata is
    // - names     ndims*CHARACTER(string_length)
    // - values    ndims*DATATYPE

    SDF_COMMON_INFO();

    SDF_READ_ENTRY_ARRAY_STRING(b->material_names, b->ndims);
    switch (b->datatype) {
    case(SDF_DATATYPE_INTEGER4):
        SDF_READ_ENTRY_ARRAY_INT4(i4, b->ndims);
        b->data = i4;
        break;
    case(SDF_DATATYPE_INTEGER8):
        SDF_READ_ENTRY_ARRAY_INT8(i8, b->ndims);
        b->data = i8;
        break;
    case(SDF_DATATYPE_REAL4):
        SDF_READ_ENTRY_ARRAY_REAL4(r4, b->ndims);
        b->data = r4;
        break;
    case(SDF_DATATYPE_REAL8):
        SDF_READ_ENTRY_ARRAY_REAL8(r8, b->ndims);
        b->data = r8;
        break;
    case(SDF_DATATYPE_LOGICAL):
        SDF_READ_ENTRY_ARRAY_LOGICAL(logical, b->ndims);
        b->data = logical;
        break;
    case(SDF_DATATYPE_CHARACTER):
        SDF_READ_ENTRY_ARRAY_STRING(string, b->ndims);
        b->data = string;
        break;
    }
    b->nmaterial_names = b->ndims;
    b->done_data = 1;

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



static int sdf_read_array(sdf_file_t *h)
{
    sdf_block_t *b = h->current_block;
    int n;

    if (b->done_data) return 0;
    if (!b->done_info) sdf_read_plain_variable_info(h);

    h->current_location = b->data_location;

    sdf_array_datatype(h);

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



static int sdf_read_datablock(sdf_file_t *h)
{
    sdf_block_t *b = h->current_block;
#ifndef PARALLEL
    size_t mlen, mstart, moff;
#endif

    if (b->done_data) return 0;
    if (!b->done_info) sdf_read_datablock_info(h);

    h->current_location = b->data_location;

#if !defined(PARALLEL) && !defined(_WIN32)
    if (h->mmap) {
        mlen = sysconf(_SC_PAGESIZE);
        mstart = mlen * (h->current_location / mlen);
        moff = h->current_location - mstart;
        b->mmap_len = mlen = b->data_length + moff;
        b->mmap = mmap(NULL, mlen, PROT_READ, MAP_SHARED, h->fd, mstart);
        b->data = moff + b->mmap;
    } else
#endif
    {
        if (b->data) free(b->data);
        b->data = malloc(b->data_length);
        sdf_seek(h);
        sdf_read_bytes(h, b->data, b->data_length);
    }

    b->done_data = 1;

    return 0;
}
