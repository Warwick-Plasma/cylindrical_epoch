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
#include <ctype.h>
#include <sdf.h>
#include "sdf_output.h"
#include "sdf_control.h"
#include "sdf_util.h"

#ifdef PARALLEL
# include <mpi.h>
#endif

#ifdef _WIN32
# include <io.h>
# define FTRUNCATE _chsize
#else
# include <unistd.h>
# define FTRUNCATE ftruncate
#endif

/*
 * @defgroup output
 * @brief Routines for writing to an SDF file
 */

static int write_constant(sdf_file_t *h);
static int write_stitched(sdf_file_t *h);
static int write_stitched_material(sdf_file_t *h);
static int write_stitched_matvar(sdf_file_t *h);
static int write_stitched_species(sdf_file_t *h);
static int write_run_info_meta(sdf_file_t *h);
static int write_data(sdf_file_t *h);



int sdf_write_bytes(sdf_file_t *h, void *buf, int buflen)
{
#ifdef PARALLEL
    return MPI_File_write(h->filehandle, buf, buflen, MPI_BYTE,
            MPI_STATUS_IGNORE);
#else
    return (1 != fwrite(buf, buflen, 1, h->filehandle));
#endif
}



int sdf_write_at(sdf_file_t *h, off_t offset, void *buf, int buflen)
{
    sdf_seek_set(h, offset);
#ifdef PARALLEL
    return MPI_File_write(h->filehandle, buf, buflen, MPI_BYTE,
            MPI_STATUS_IGNORE);
#else
    return (1 != fwrite(buf, buflen, 1, h->filehandle));
#endif
}



int sdf_flush(sdf_file_t *h)
{
    //errcode += sdf_update(h);
#ifdef PARALLEL
    return MPI_File_sync(h->filehandle);
#else
    return fflush(h->filehandle);
#endif
}



static int sdf_truncate(sdf_file_t *h, int64_t size)
{
#ifdef PARALLEL
    MPI_Offset length = size;
    return MPI_File_set_size(h->filehandle, length);
#else
    off_t length = size;
    return FTRUNCATE(fileno(h->filehandle), length);
#endif
}



static size_t trimwhitespace(const char *str_in, char *str_out, size_t len)
{
    size_t out_size, i1, i2;
    const char *end;

    if (len == 0) return 0;

    // Trim leading space
    while (isspace(*str_in)) str_in++;

    if (*str_in == 0) { // All spaces?
        *str_out = 0;
        return 1;
    }

    // Trim trailing space
    end = str_in + strlen(str_in) - 1;
    while (end > str_in && isspace(*end)) end--;
    end++;

    // Set output size to minimum of trimmed string length and
    // buffer size minus 1
    i1 = end - str_in;
    i2 = len - 1;
    out_size = i1 < i2 ? i1 : i2;

    // Copy trimmed string and add null terminator
    memcpy(str_out, str_in, out_size);
    str_out[out_size] = 0;

    return out_size;
}



static int safe_copy_string(char *s1, char *s2)
{
    int len1, len2;

    len1 = strlen(s1);
    len2 = strlen(s2);

    memset(s2, 0, len2);
    memcpy(s2, s1, len1);

    return 0;
}



static int sdf_safe_write_string_len(sdf_file_t *h, char *string, int length)
{
    char *output;
    int len_s, errcode = 0;

    output = malloc(length);
    len_s = trimwhitespace(string, output, length);

    if (len_s > length && h->rank == h->rank_master) {
        printf("*** WARNING ***\n");
        printf("Output string \"%s\" has been truncated.", output);
    }

    if (len_s+1 < length) output[len_s+1] = 0;

    errcode = sdf_write_bytes(h, output, length);
    free(output);

    return errcode;
}



static int sdf_safe_write_string(sdf_file_t *h, char *string)
{
    return sdf_safe_write_string_len(h, string, h->string_length);
}



static int sdf_safe_write_id(sdf_file_t *h, char *string)
{
    int length = h->id_length;
    return sdf_safe_write_string_len(h, string, length);
}



static int write_block(sdf_file_t *h)
{
    sdf_block_t *b;

    b = h->current_block;
    b->block_start = h->current_location;

    switch (b->blocktype) {
    case SDF_BLOCKTYPE_CONSTANT:
        write_constant(h);
        break;
    case SDF_BLOCKTYPE_STITCHED_MATERIAL:
    case SDF_BLOCKTYPE_CONTIGUOUS_MATERIAL:
        write_stitched_material(h);
        break;
    case SDF_BLOCKTYPE_STITCHED_MATVAR:
    case SDF_BLOCKTYPE_CONTIGUOUS_MATVAR:
        write_stitched_matvar(h);
        break;
    case SDF_BLOCKTYPE_STITCHED:
    case SDF_BLOCKTYPE_CONTIGUOUS:
    case SDF_BLOCKTYPE_STITCHED_TENSOR:
    case SDF_BLOCKTYPE_CONTIGUOUS_TENSOR:
        write_stitched(h);
        break;
    case SDF_BLOCKTYPE_STITCHED_SPECIES:
    case SDF_BLOCKTYPE_CONTIGUOUS_SPECIES:
        write_stitched_species(h);
        break;
    case SDF_BLOCKTYPE_RUN_INFO:
        write_run_info_meta(h);
        break;
    case SDF_BLOCKTYPE_ARRAY:
    case SDF_BLOCKTYPE_PLAIN_MESH:
    case SDF_BLOCKTYPE_PLAIN_VARIABLE:
    case SDF_BLOCKTYPE_CPU_SPLIT:
        write_data(h);
        break;
    default:
        printf("WARNING! Ignored id: %s\n", b->id);
    }
#if 0
    SDF_BLOCKTYPE_POINT_MESH,
    SDF_BLOCKTYPE_POINT_VARIABLE,
    SDF_BLOCKTYPE_SOURCE,
    SDF_BLOCKTYPE_SPECIES,
    SDF_BLOCKTYPE_PLAIN_DERIVED,
    SDF_BLOCKTYPE_POINT_DERIVED,
    SDF_BLOCKTYPE_STITCHED_OBSTACLE_GROUP,
    SDF_BLOCKTYPE_UNSTRUCTURED_MESH,
    SDF_BLOCKTYPE_LAGRANGIAN_MESH,
    SDF_BLOCKTYPE_STATION,
    SDF_BLOCKTYPE_STATION_DERIVED,
#endif
    return 0;
}



static int write_header(sdf_file_t *h)
{
    int32_t int4;
    int errcode;
    char padding[5];

    errcode = 0;

    if (h->done_header) {
        if (h->rank == h->rank_master) {
            printf("*** WARNING ***\n");
            printf("SDF header already written. Ignoring extra call.\n");
        }
        return errcode;
    }

    // Currently no blocks written
    h->nblocks = 0;
    h->summary_location = h->first_block_location;
    h->summary_size = 0;
    h->current_location = 0;
    //h->data_location = 0;

    if (h->rank == h->rank_master) {
        errcode += sdf_seek(h);

        // Write the header
        errcode += sdf_write_bytes(h, SDF_MAGIC, 4);

        int4 = SDF_ENDIANNESS;
        errcode += sdf_write_bytes(h, &int4, SOI4);

        int4 = SDF_VERSION;
        errcode += sdf_write_bytes(h, &int4, SOI4);

        int4 = SDF_REVISION;
        errcode += sdf_write_bytes(h, &int4, SOI4);

        errcode += sdf_safe_write_id(h, h->code_name);

        errcode += sdf_write_bytes(h, &h->first_block_location, SOI8);

        // Must be consistent with the c_summary_offset in sdf_common
        errcode += sdf_write_bytes(h, &h->summary_location, SOI8);

        errcode += sdf_write_bytes(h, &h->summary_size, SOI4);

        errcode += sdf_write_bytes(h, &h->nblocks, SOI4);

        errcode += sdf_write_bytes(h, &h->block_header_length, SOI4);

        errcode += sdf_write_bytes(h, &h->step, SOI4);

        errcode += sdf_write_bytes(h, &h->time, SOF8);

        errcode += sdf_write_bytes(h, &h->jobid1, SOI4);

        errcode += sdf_write_bytes(h, &h->jobid2, SOI4);

        errcode += sdf_write_bytes(h, &h->string_length, SOI4);

        errcode += sdf_write_bytes(h, &h->code_io_version, SOI4);

        errcode += sdf_write_bytes(h, &h->restart_flag, 1);

        errcode += sdf_write_bytes(h, &h->other_domains, 1);

        errcode += sdf_write_bytes(h, &h->station_file, 1);

        memset(padding, 0, 5);
        errcode += sdf_write_bytes(h, padding, 5);
    }

    h->current_location = h->first_block_location;
    h->done_header = 1;

    return 0;
}



static int write_block_header(sdf_file_t *h)
{
    int errcode = 0, block_info_length;
    sdf_block_t *b = h->current_block;

    if (!b || b->done_header || !b->in_file) return errcode;

    if (h->use_summary) {
        h->current_location = b->block_start = b->summary_block_start;
        b->next_block_location = b->summary_next_block_location;
    } else {
        h->current_location = b->block_start = b->inline_block_start;
        b->next_block_location = b->inline_next_block_location;
    }

    // If this routine is changed then the value of h->block_header_length
    // must be changed accordingly in sdf_write_header

    // Write the block header
    if (h->rank == h->rank_master) {
        errcode += sdf_seek(h);

        errcode += sdf_write_bytes(h, &b->next_block_location, SOI8);

        errcode += sdf_write_bytes(h, &b->data_location, SOI8);

        errcode += sdf_safe_write_id(h, b->id);

        errcode += sdf_write_bytes(h, &b->data_length, SOI8);

        errcode += sdf_write_bytes(h, &b->blocktype, SOI4);

        errcode += sdf_write_bytes(h, &b->datatype, SOI4);

        errcode += sdf_write_bytes(h, &b->ndims, SOI4);

        errcode += sdf_safe_write_string(h, b->name);

        block_info_length = b->info_length - h->block_header_length;
        errcode += sdf_write_bytes(h, &block_info_length, SOI4);
    }

    h->current_location = b->block_start + h->block_header_length;
    b->done_header = 1;

    return errcode;
}



static int write_constant(sdf_file_t *h)
{
    int errcode;
    sdf_block_t *b = h->current_block;

    if (!b || !b->in_file) return 0;

    b->nelements = b->ndims = 1;

    // Metadata is
    // - value     SDF_TYPE_SIZES[b->datatype]

    b->info_length = h->block_header_length + SDF_TYPE_SIZES[b->datatype];
    b->data_length = 0;

    // Write header
    errcode = write_block_header(h);

    // Write data (in metadata section)
    if (h->rank == h->rank_master) {
        errcode += sdf_write_bytes(h, &b->const_value,
                SDF_TYPE_SIZES[b->datatype]);
    }

    h->current_location = b->block_start + b->info_length;
    b->done_info = 1;
    b->done_data = 1;

    return errcode;
}



static int write_array_meta(sdf_file_t *h)
{
    int errcode, i;
    int32_t int4;
    sdf_block_t *b = h->current_block;

    if (!b || !b->in_file) return 0;

    b->nelements = 1;
    for (i=0; i < b->ndims; i++)
        b->nelements *= b->dims[i];

    // Metadata is
    // - dims      ndims*INTEGER(i4)

    b->info_length = h->block_header_length + b->ndims * SOI4;
    b->data_length = b->nelements * SDF_TYPE_SIZES[b->datatype];

    // Write header
    errcode = write_block_header(h);

    // Write metadata
    if (h->rank == h->rank_master) {
        for (i=0; i < b->ndims; i++) {
            int4 = b->dims[i];
            errcode += sdf_write_bytes(h, &int4, SOI4);
        }
    }

    h->current_location = b->block_start + b->info_length;
    b->done_info = 1;

    return errcode;
}



static int write_cpu_split_meta(sdf_file_t *h)
{
    int errcode, i;
    int32_t int4;
    sdf_block_t *b = h->current_block;

    if (!b || !b->in_file) return 0;

    b->blocktype = SDF_BLOCKTYPE_CPU_SPLIT;

    // Metadata is
    // - type      INTEGER(i4)
    // - dims      ndims*INTEGER(i4)

    b->info_length = h->block_header_length + (b->ndims + 1) * SOI4;
    b->data_length = b->nelements * SDF_TYPE_SIZES[b->datatype];

    // Write header
    errcode = write_block_header(h);

    // Write metadata
    if (h->rank == h->rank_master) {
        errcode += sdf_write_bytes(h, &b->geometry, SOI4);

        for (i=0; i < b->ndims; i++) {
            int4 = b->dims[i];
            errcode += sdf_write_bytes(h, &int4, SOI4);
        }
    }

    h->current_location = b->block_start + b->info_length;
    b->done_info = 1;

    return errcode;
}



static int write_stitched(sdf_file_t *h)
{
    int errcode, i;
    sdf_block_t *b = h->current_block;

    if (!b || !b->in_file) return 0;

    // Metadata is
    // - stagger   INTEGER(i4)
    // - meshid    CHARACTER(id_length)
    // - varids    ndims*CHARACTER(id_length)

    b->info_length = h->block_header_length + SOI4
            + (b->ndims + 1) * h->id_length;
    if (b->blocktype == SDF_BLOCKTYPE_STITCHED
            || b->blocktype == SDF_BLOCKTYPE_STITCHED_TENSOR)
        b->data_length = 0;

    // Write header
    errcode = write_block_header(h);

    // Write metadata
    if (h->rank == h->rank_master) {
        errcode += sdf_write_bytes(h, &b->stagger, SOI4);

        errcode += sdf_safe_write_id(h, b->mesh_id);

        for (i=0; i < b->ndims; i++)
            errcode += sdf_safe_write_id(h, b->variable_ids[i]);
    }

    h->current_location = b->block_start + b->info_length;

    b->done_info = 1;
    b->done_data = 1;

    return errcode;
}



static int write_stitched_material(sdf_file_t *h)
{
    int errcode, i;
    sdf_block_t *b = h->current_block;

    if (!b || !b->in_file) return 0;

    // Metadata is
    // - stagger   INTEGER(i4)
    // - meshid    CHARACTER(id_length)
    // - material_names ndims*CHARACTER(string_length)
    // - varids    ndims*CHARACTER(id_length)

    //b->info_length = h->block_header_length + SOI4
    //        + (b->ndims + 1) * h->id_length + b->ndims * h->string_length;
    if (b->blocktype == SDF_BLOCKTYPE_STITCHED_MATERIAL) b->data_length = 0;

    // Write header
    errcode = write_block_header(h);

    // Write metadata
    if (h->rank == h->rank_master) {
        errcode += sdf_write_bytes(h, &b->stagger, SOI4);

        errcode += sdf_safe_write_id(h, b->mesh_id);

        for (i=0; i < b->ndims; i++)
            errcode += sdf_safe_write_string(h, b->material_names[i]);

        for (i=0; i < b->ndims; i++)
            errcode += sdf_safe_write_id(h, b->variable_ids[i]);
    }

/*
    if (b->data_length > 0)
        h->current_location = b->next_block_location;
    else
        h->current_location = b->block_start + b->info_length;
*/

    b->done_info = 1;
    b->done_data = 1;

    return errcode;
}



static int write_stitched_matvar(sdf_file_t *h)
{
    int errcode, i;
    sdf_block_t *b = h->current_block;

    if (!b || !b->in_file) return 0;

    // Metadata is
    // - stagger   INTEGER(i4)
    // - meshid    CHARACTER(id_length)
    // - matid     CHARACTER(id_length)
    // - varids    ndims*CHARACTER(id_length)

    b->info_length = h->block_header_length + SOI4
            + (b->ndims + 2) * h->id_length;
    if (b->blocktype == SDF_BLOCKTYPE_STITCHED_MATVAR) b->data_length = 0;

    // Write header
    errcode = write_block_header(h);

    // Write metadata
    if (h->rank == h->rank_master) {
        errcode += sdf_write_bytes(h, &b->stagger, SOI4);

        errcode += sdf_safe_write_id(h, b->mesh_id);

        errcode += sdf_safe_write_id(h, b->material_id);

        for (i=0; i < b->ndims; i++)
            errcode += sdf_safe_write_id(h, b->variable_ids[i]);
    }

/*
    if (b->data_length > 0)
        h->current_location = b->next_block_location;
    else
        h->current_location = b->block_start + b->info_length;
*/

    h->current_location = b->block_start + b->info_length;
    b->done_info = 1;
    b->done_data = 1;

    return errcode;
}



static int write_stitched_species(sdf_file_t *h)
{
    int errcode, i;
    sdf_block_t *b = h->current_block;

    if (!b || !b->in_file) return 0;

    // Metadata is
    // - stagger   INTEGER(i4)
    // - meshid    CHARACTER(id_length)
    // - matid     CHARACTER(id_length)
    // - matname   CHARACTER(string_length)
    // - specnames ndims*CHARACTER(string_length)
    // - varids    ndims*CHARACTER(id_length)

    b->info_length = h->block_header_length + SOI4
            + (b->ndims + 2) * h->id_length + (b->ndims + 1) * h->string_length;
    if (b->blocktype == SDF_BLOCKTYPE_STITCHED_SPECIES) b->data_length = 0;

    // Write header
    errcode = write_block_header(h);

    // Write metadata
    if (h->rank == h->rank_master) {
        errcode += sdf_write_bytes(h, &b->stagger, SOI4);

        errcode += sdf_safe_write_id(h, b->mesh_id);

        errcode += sdf_safe_write_id(h, b->material_id);

        errcode += sdf_safe_write_string(h, b->material_name);

        for (i=0; i < b->ndims; i++)
            errcode += sdf_safe_write_string(h, b->material_names[i]);

        for (i=0; i < b->ndims; i++)
            errcode += sdf_safe_write_id(h, b->variable_ids[i]);
    }

    if (b->data_length > 0)
        h->current_location = b->next_block_location;
    else
        h->current_location = b->block_start + b->info_length;

    b->done_info = 1;
    b->done_data = 1;

    return errcode;
}



static int write_stitched_obstacle_group(sdf_file_t *h)
{
    int errcode, i;
    sdf_block_t *b = h->current_block;

    if (!b || !b->in_file) return 0;

    // Metadata is
    // - stagger        INTEGER(i4)
    // - obstacle_id    CHARACTER(id_length)
    // - vfm_id         CHARACTER(id_length)
    // - obstacle_names ndims*CHARACTER(string_length)

    b->info_length = h->block_header_length + SOI4 + 2 * h->id_length
            + b->ndims * h->string_length;
    b->data_length = 0;

    // Write header
    errcode = write_block_header(h);

    // Write metadata
    if (h->rank == h->rank_master) {
        errcode += sdf_write_bytes(h, &b->stagger, SOI4);

        errcode += sdf_safe_write_id(h, b->obstacle_id);

        errcode += sdf_safe_write_id(h, b->vfm_id);

        for (i=0; i < b->ndims; i++)
            errcode += sdf_safe_write_string(h, b->material_names[i]);
    }

    h->current_location = b->block_start + b->info_length;
    b->done_info = 1;

    return errcode;
}



static int write_run_info_meta(sdf_file_t *h)
{
    int errcode;
    sdf_block_t *b = h->current_block;
    struct run_info *run = b->data;

    if (!b || !b->in_file) return 0;

    b->datatype = SDF_DATATYPE_OTHER;
    b->blocktype = SDF_BLOCKTYPE_RUN_INFO;
    b->ndims = 1;

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

    b->info_length = h->block_header_length + 6 * SOI4 + SOI8
            + 4 * h->string_length;
    b->data_length = 0;

    // Write header
    errcode = write_block_header(h);

    // Write metadata
    if (h->rank == h->rank_master) {
        errcode += sdf_write_bytes(h, &run->version, SOI4);

        errcode += sdf_write_bytes(h, &run->revision, SOI4);

        errcode += sdf_safe_write_string(h, run->commit_id);

        errcode += sdf_safe_write_string(h, run->sha1sum);

        errcode += sdf_safe_write_string(h, run->compile_machine);

        errcode += sdf_safe_write_string(h, run->compile_flags);

        errcode += sdf_write_bytes(h, &run->defines, SOI8);

        errcode += sdf_write_bytes(h, &run->compile_date, SOI4);

        errcode += sdf_write_bytes(h, &run->run_date, SOI4);

        errcode += sdf_write_bytes(h, &run->io_date, SOI4);

        errcode += sdf_write_bytes(h, &run->minor_rev, SOI4);
    }

    h->current_location = b->block_start + b->info_length;
    b->done_info = 1;
    b->done_data = 1;

    return errcode;
}



static int write_plain_mesh_meta(sdf_file_t *h)
{
    int errcode, i;
    int32_t int4;
    sdf_block_t *b = h->current_block;

    if (!b || !b->in_file) return 0;

/*
    if (b->blocktype == SDF_BLOCKTYPE_LAGRANGIAN_MESH) {
        b->nelements = b->ndims;
        for (i=0; i < b->ndims; i++)
            b->nelements *= b->dims[i];
    } else {
        b->blocktype = SDF_BLOCKTYPE_PLAIN_MESH;
        b->nelements = 0;
        for (i=0; i < b->ndims; i++)
            b->nelements += b->dims[i];
    }
*/

    // Metadata is
    // - mults     REAL(r8), DIMENSION(ndims)
    // - labels    CHARACTER(id_length), DIMENSION(ndims)
    // - units     CHARACTER(id_length), DIMENSION(ndims)
    // - geometry  INTEGER(i4)
    // - minval    REAL(r8), DIMENSION(ndims)
    // - maxval    REAL(r8), DIMENSION(ndims)
    // - dims      INTEGER(i4), DIMENSION(ndims)

    b->info_length = h->block_header_length + (b->ndims + 1) * SOI4
            + (3 * b->ndims) * SOF8 + 2 * b->ndims * h->id_length;
    b->data_length = b->nelements * SDF_TYPE_SIZES[b->datatype];

    // Write header
    errcode = write_block_header(h);

    // Write metadata
    if (h->rank == h->rank_master) {
        errcode += sdf_write_bytes(h, b->dim_mults, b->ndims * SOF8);

        for (i=0; i < b->ndims; i++)
            errcode += sdf_safe_write_id(h, b->dim_labels[i]);

        for (i=0; i < b->ndims; i++)
            errcode += sdf_safe_write_id(h, b->dim_units[i]);

        errcode += sdf_write_bytes(h, &b->geometry, SOI4);

        errcode += sdf_write_bytes(h, b->extents, 2 * b->ndims * SOF8);

        for (i=0; i < b->ndims; i++) {
            int4 = b->dims[i];
            errcode += sdf_write_bytes(h, &int4, SOI4);
        }
    }

    h->current_location = b->block_start + b->info_length;
    b->done_info = 1;

    return errcode;
}



static int write_plain_variable_meta(sdf_file_t *h)
{
    int errcode, i;
    int32_t int4;
    sdf_block_t *b = h->current_block;

    if (!b || !b->in_file) return 0;

    b->blocktype = SDF_BLOCKTYPE_PLAIN_VARIABLE;

    b->nelements = 1;
    for (i=0; i < b->ndims; i++)
        b->nelements *= b->dims[i];

    // Metadata is
    // - mult      REAL(r8)
    // - units     CHARACTER(id_length)
    // - meshid    CHARACTER(id_length)
    // - dims      INTEGER(i4), DIMENSION(ndims)
    // - stagger   INTEGER(i4)

    b->info_length = h->block_header_length + (b->ndims + 1) * SOI4 + SOF8
            + 2 * h->id_length;
    b->data_length = b->nelements * SDF_TYPE_SIZES[b->datatype];

    // Write header
    errcode = write_block_header(h);

    // Write metadata
    if (h->rank == h->rank_master) {
        errcode += sdf_seek(h);

        errcode += sdf_write_bytes(h, &b->mult, SOF8);

        errcode += sdf_safe_write_id(h, b->units);

        errcode += sdf_safe_write_id(h, b->mesh_id);

        for (i=0; i < b->ndims; i++) {
            int4 = b->dims[i];
            errcode += sdf_write_bytes(h, &int4, SOI4);
        }

        errcode += sdf_write_bytes(h, &b->stagger, SOI4);
    }

    h->current_location = b->block_start + b->info_length;
    b->done_info = 1;

    return errcode;
}



static int write_point_mesh_meta(sdf_file_t *h)
{
    int errcode, i;
    sdf_block_t *b = h->current_block;

    if (!b || !b->in_file) return 0;

    //b->blocktype = SDF_BLOCKTYPE_POINT_MESH;

    // Metadata is
    // - mults     REAL(r8), DIMENSION(ndims)
    // - labels    CHARACTER(id_length), DIMENSION(ndims)
    // - units     CHARACTER(id_length), DIMENSION(ndims)
    // - geometry  INTEGER(i4)
    // - minval    REAL(r8), DIMENSION(ndims)
    // - maxval    REAL(r8), DIMENSION(ndims)
    // - npoints   INTEGER(i8)
    // - speciesid CHARACTER(id_length)

    b->info_length = h->block_header_length + SOI4 + SOI8
            + (3 * b->ndims) * SOF8 + (2 * b->ndims + 1) * h->id_length;
    b->data_length = b->nelements * SDF_TYPE_SIZES[b->datatype];

    // Write header
    errcode = write_block_header(h);

    // Write metadata
    if (h->rank == h->rank_master) {
        errcode += sdf_write_bytes(h, b->dim_mults, b->ndims * SOF8);

        for (i=0; i < b->ndims; i++)
            errcode += sdf_safe_write_id(h, b->dim_labels[i]);

        for (i=0; i < b->ndims; i++)
            errcode += sdf_safe_write_id(h, b->dim_units[i]);

        errcode += sdf_write_bytes(h, &b->geometry, SOI4);

        errcode += sdf_write_bytes(h, b->extents, 2 * b->ndims * SOF8);

        errcode += sdf_write_bytes(h, &b->nelements, SOI8);

        errcode += sdf_safe_write_id(h, b->material_id);
    }

    h->current_location = b->block_start + b->info_length;
    b->done_info = 1;

    return errcode;
}



static int write_point_variable_meta(sdf_file_t *h)
{
    int errcode;
    sdf_block_t *b = h->current_block;

    if (!b || !b->in_file) return 0;

    b->blocktype = SDF_BLOCKTYPE_POINT_VARIABLE;

    // Metadata is
    // - mult      REAL(r8)
    // - units     CHARACTER(id_length)
    // - meshid    CHARACTER(id_length)
    // - npoints   INTEGER(i8)
    // - speciesid CHARACTER(id_length)

    b->info_length = h->block_header_length + SOI8 + SOF8 + 3 * h->id_length;
    b->data_length = b->nelements * SDF_TYPE_SIZES[b->datatype];

    // Write header
    errcode = write_block_header(h);

    // Write metadata
    if (h->rank == h->rank_master) {
        errcode += sdf_write_bytes(h, &b->mult, SOF8);

        errcode += sdf_safe_write_id(h, b->units);

        errcode += sdf_safe_write_id(h, b->mesh_id);

        errcode += sdf_write_bytes(h, &b->nelements, SOI8);

        errcode += sdf_safe_write_id(h, b->material_id);
    }

    h->current_location = b->block_start + b->info_length;
    b->done_info = 1;

    return errcode;
}



static int write_meta(sdf_file_t *h)
{
    int return_value = 0;
    sdf_block_t *b = h->current_block;

    if (!b || !b->in_file) return 0;

    b->info_length += h->block_header_length;

    switch (b->blocktype) {
    case SDF_BLOCKTYPE_PLAIN_MESH:
    case SDF_BLOCKTYPE_LAGRANGIAN_MESH:
        return_value = write_plain_mesh_meta(h);
        break;
    case SDF_BLOCKTYPE_POINT_MESH:
        return_value = write_point_mesh_meta(h);
        break;
    case SDF_BLOCKTYPE_PLAIN_VARIABLE:
        return_value = write_plain_variable_meta(h);
        break;
    case SDF_BLOCKTYPE_POINT_VARIABLE:
        return_value = write_point_variable_meta(h);
        break;
    case SDF_BLOCKTYPE_CONSTANT:
        return_value = write_constant(h);
        break;
    case SDF_BLOCKTYPE_ARRAY:
        return_value = write_array_meta(h);
        break;
    case SDF_BLOCKTYPE_CPU_SPLIT:
        return_value = write_cpu_split_meta(h);
        break;
    case SDF_BLOCKTYPE_RUN_INFO:
        return_value = write_run_info_meta(h);
        break;
    case SDF_BLOCKTYPE_SOURCE:
        return_value = write_block_header(h);
        break;
    case SDF_BLOCKTYPE_STITCHED:
    case SDF_BLOCKTYPE_CONTIGUOUS:
    case SDF_BLOCKTYPE_STITCHED_TENSOR:
    case SDF_BLOCKTYPE_CONTIGUOUS_TENSOR:
        return_value = write_stitched(h);
        break;
    case SDF_BLOCKTYPE_STITCHED_MATERIAL:
    case SDF_BLOCKTYPE_CONTIGUOUS_MATERIAL:
        return_value = write_stitched_material(h);
        break;
    case SDF_BLOCKTYPE_STITCHED_MATVAR:
    case SDF_BLOCKTYPE_CONTIGUOUS_MATVAR:
        return_value = write_stitched_matvar(h);
        break;
    case SDF_BLOCKTYPE_STITCHED_SPECIES:
    case SDF_BLOCKTYPE_CONTIGUOUS_SPECIES:
        return_value = write_stitched_species(h);
        break;
    case SDF_BLOCKTYPE_STITCHED_OBSTACLE_GROUP:
        return_value = write_stitched_obstacle_group(h);
        break;
    default:
        return_value = write_block_header(h);
        break;
    }

    h->current_location = b->block_start + b->info_length;
    b->info_length -= h->block_header_length;

    return return_value;
}



int sdf_write_meta(sdf_file_t *h)
{
    return write_meta(h);
}



static int write_data(sdf_file_t *h)
{
    int errcode = 0, i;
    sdf_block_t *b = h->current_block;

    if (!b || !b->in_file) return 0;

    // Write header
    errcode = write_meta(h);

    if (h->rank == h->rank_master) {
        errcode += sdf_seek_set(h, b->data_location);

        // Actual array
        if (b->data) {
            errcode += sdf_write_bytes(h, b->data, b->data_length);
        } else if (b->grids) {
            for (i=0; i < b->ngrids; i++)
                errcode += sdf_write_bytes(h, b->grids[i],
                        b->dims[i] * SDF_TYPE_SIZES[b->datatype]);
        }
    }

    h->current_location = b->data_location + b->data_length;
    b->done_data = 1;

    return errcode;
}



int64_t sdf_write_new_summary(sdf_file_t *h)
{
    sdf_block_t *b, *next;
    int64_t total_summary_size = 0, extent = h->first_block_location;
    int64_t summary_size, sz;
    int errcode = 0, use_summary;

    // First find the furthest extent of the inline metadata and/or data
    next = h->blocklist;
    while (next) {
        b = next;
        next = next->next;
        if (!b->in_file) continue;

        sz = b->inline_block_start + h->block_header_length + b->info_length;
        if (sz > extent) extent = sz;

        sz = b->data_location + b->data_length;
        if (sz > extent) extent = sz;
    }

    h->current_location = h->summary_location = extent;

    use_summary = h->use_summary;
    h->use_summary = 1;
    next = h->blocklist;
    while (next) {
        b = h->current_block = next;
        next = next->next;
        if (!b->in_file) continue;
        b->done_header = 0;
        summary_size = h->block_header_length + b->info_length;
        total_summary_size += summary_size;
        b->summary_block_start = b->block_start = h->current_location;
        b->summary_next_block_location = b->next_block_location
                = b->block_start + summary_size;
        errcode += write_meta(h);
    }
    sdf_truncate(h, h->current_location);
    h->use_summary = use_summary;

    h->summary_size = total_summary_size;

    return total_summary_size;
}



static int sdf_write_header(sdf_file_t *h, char *code_name, int code_io_version,
        int step, double time, char restart, char station_file, int jobid1,
        int jobid2)
{
    int errcode = 0;

    if (h->code_name) free(h->code_name);
    errcode += safe_copy_string(code_name, h->code_name);

    h->step = step;
    h->time = time;
    h->restart_flag = restart;
    h->jobid1 = jobid1;
    h->jobid2 = jobid2;
    h->station_file = station_file;

    return write_header(h);
}



int sdf_write(sdf_file_t *h)
{
    int errcode;
    sdf_block_t *b;

    errcode = write_header(h);
    sdf_flush(h);

    b = h->blocklist;
    while (b) {
        h->current_block = b;
        b->done_header = b->done_info = b->done_data = 0;
        errcode += write_block(h);
        //sdf_flush(h);
        b = b->next;
    }

    return errcode;
}
