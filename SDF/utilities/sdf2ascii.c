/*
 * sdf2ascii - print the contents of an SDF file
 * Copyright (C) 2010-2016 SDF Development Team
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <sys/stat.h>
#include "sdf.h"
#include "commit_info.h"

#define SDF_DEBUG

#ifdef PARALLEL
#include <mpi.h>
#endif

#define VERSION "2.6.7"

#define DBG_FLUSH() do { \
        if (h && h->dbg_buf) { \
            h->dbg = h->dbg_buf; *h->dbg = '\0'; \
        } \
    } while (0)

#define DBG_PRINT_FLUSH() do { \
        if (h && h->dbg_buf) { \
            printf("%s", h->dbg_buf); h->dbg = h->dbg_buf; *h->dbg = '\0'; \
        } \
    } while (0)

int metadata, contents, debug, single, use_mmap, ignore_summary;
int element_count, ignore_nblocks;
struct id_list {
    char *id;
    struct id_list *next;
} *variable_ids, *last_id;

int nrange;

struct range_type {
    int start, end;
} *range_list;


void usage(int err)
{
    fprintf(stderr, "usage: sdf2ascii [options] <sdf_filename>\n");
    fprintf(stderr, "\noptions:\n\
  -h --help            Show this usage message\n\
  -n --no-metadata     Don't show metadata blocks (shown by default)\n\
  -c --contents        Show block's data content\n\
  -s --single          Convert block data to single precision\n\
  -v --variable=id     Find the block with id matching 'id'\n\
  -m --mmap            Use mmap'ed file I/O\n\
  -i --no-summary      Ignore the metadata summary\n\
  -b --no-nblocks      Ignore the header value for nblocks\n\
  -C --count=n         When printing array contents, write 'n' elements per\n\
                       line.\n\
  -V --version         Print version information and exit\n\
");
/*
  -D --debug           Show the contents of the debug buffer\n\
*/

    exit(err);
}


int range_sort(const void *v1, const void *v2)
{
    struct range_type *a = (struct range_type *)v1;
    struct range_type *b = (struct range_type *)v2;

    return (a->start - b->start);
}


char *parse_args(int *argc, char ***argv)
{
    char *ptr, *file = NULL;
    int c, i, err, range, sz, nbuf;
    struct range_type *range_tmp;
    struct stat statbuf;
    static struct option longopts[] = {
        { "no-nblocks",    no_argument,       NULL, 'b' },
        { "contents",      no_argument,       NULL, 'c' },
        { "count",         required_argument, NULL, 'C' },
        { "help",          no_argument,       NULL, 'h' },
        { "no-summary",    no_argument,       NULL, 'i' },
        { "mmap",          no_argument,       NULL, 'm' },
        { "no-metadata",   no_argument,       NULL, 'n' },
        { "single",        no_argument,       NULL, 's' },
        { "variable",      required_argument, NULL, 'v' },
        { "version",       no_argument,       NULL, 'V' },
        { NULL,            0,                 NULL,  0  }
        //{ "debug",         no_argument,       NULL, 'D' },
    };

    metadata = debug = 1;
    contents = single = use_mmap = ignore_summary = ignore_nblocks = 0;
    variable_ids = NULL;
    last_id = NULL;
    nbuf = nrange = element_count = 0;
    sz = sizeof(struct range_type);

    while ((c = getopt_long(*argc, *argv,
            "bcC:himnsv:V", longopts, NULL)) != -1) {
        switch (c) {
        case 'b':
            ignore_nblocks = 1;
            break;
        case 'c':
            contents = 1;
            break;
        case 'C':
            element_count = strtol(optarg, NULL, 10);
            if (element_count < 1) element_count = 1;
            break;
        case 'h':
            usage(0);
            break;
        case 'i':
            ignore_summary = 1;
            break;
        case 'm':
            use_mmap = 1;
            break;
        case 'n':
            metadata = 0;
            break;
        case 's':
            single = 1;
            break;
        case 'V':
            printf("sdf2ascii version %s\n", VERSION);
            printf("commit info: %s, %s\n", SDF_COMMIT_ID, SDF_COMMIT_DATE);
            printf("library commit info: %s, %s\n",
                   sdf_get_library_commit_id(), sdf_get_library_commit_date());
            exit(0);
            break;
        case 'v':
            if (*optarg >= '0' && *optarg <= '9') {
                ptr = optarg;
                range = 0;
                while (ptr < optarg + strlen(optarg) + 1) {
                    if (range) {
                        i = (int)strtol(ptr, &ptr, 10);
                        if (i == 0)
                            range_list[nrange-1].end = -1;
                        else if (i < range_list[nrange-1].start)
                            nrange--;
                        else
                            range_list[nrange-1].end = i;
                        range = 0;
                    } else {
                        nrange++;
                        // Grow array if necessary
                        if (nrange > nbuf) {
                            if (nbuf == 0) {
                                nbuf = 128;
                                range_list = malloc(nbuf * sz);
                            } else {
                                i = 2 * nbuf;

                                range_tmp = malloc(i * sz);
                                memcpy(range_tmp, range_list, nbuf * sz);
                                free(range_list);
                                range_list = range_tmp;

                                nbuf = i;
                            }
                        }

                        i = (int)strtol(ptr, &ptr, 10);
                        range_list[nrange-1].start = i;
                        range_list[nrange-1].end = i;
                        if (*ptr == '-') range = 1;
                    }

                    ptr++;
                }
            } else {
                if (!variable_ids) {
                    last_id = variable_ids = malloc(sizeof(*variable_ids));
                } else {
                    last_id->next = malloc(sizeof(*variable_ids));
                    last_id = last_id->next;
                }
                last_id->next = NULL;
                last_id->id = malloc(strlen(optarg)+1);
                memcpy(last_id->id, optarg, strlen(optarg)+1);
            }
            break;
        default:
            usage(1);
        }
    }

    if ((optind+1) == *argc) {
        file = (*argv)[optind];
        err = lstat(file, &statbuf);
        if (err) {
            fprintf(stderr, "Error opening file %s\n", file);
            exit(1);
        }
    } else {
        fprintf(stderr, "No file specified\n");
        usage(1);
    }

    if (nrange > 0) {
        // Sanitize range list
        qsort(range_list, nrange, sz, &range_sort);
        for (i=1; i < nrange; ) {
            if (range_list[i].start <= range_list[i-1].end+1) {
                if (range_list[i].end > range_list[i-1].end)
                    range_list[i-1].end = range_list[i].end;
                memcpy(range_list+i, range_list+i+1, (nrange-i) * sz);
                nrange--;
            } else
                i++;
        }

        // Shrink array
        range_tmp = malloc(nrange * sz);
        memcpy(range_tmp, range_list, nrange * sz);
        free(range_list);
        range_list = range_tmp;
    }

    return file;
}


int main(int argc, char **argv)
{
    char *file = NULL;
    int i, n, block, err, found, print_block, range_start;
    sdf_file_t *h;
    sdf_block_t *b;
    comm_t comm;

    if (sdf_has_debug_info() == 0) {
        fprintf(stderr, "ERROR: sdf2ascii has been linked against a copy of "
                "the SDF C library that has\n"
                "       no debug info. Both library and sdf2ascii utility must "
                "be recompiled.\n");
        return 1;
    }

    file = parse_args(&argc, &argv);

#ifdef PARALLEL
    MPI_Init(&argc, &argv);
    MPI_Comm_dup(MPI_COMM_WORLD, &comm);
#else
    comm = 0;
#endif

    h = sdf_open(file, comm, SDF_READ, use_mmap);
    if (!h) {
        fprintf(stderr, "Error opening file %s\n", file);
        return 1;
    }
    if (element_count > 0) h->array_count = element_count;

    h->use_float = single;
    h->print = debug;
    if (ignore_summary) h->use_summary = 0;
    if (ignore_nblocks) h->ignore_nblocks = 1;

    sdf_read_header(h);
    h->current_block = NULL;

    // If nblocks is negative then the file is corrupt
    if (h->nblocks < 0) {
        block = (-h->nblocks) / 64;
        err = -h->nblocks - 64 * block;
        fprintf(stderr, "Error code %s found at block %i\n",
                sdf_error_codes_c[err], block);
        //return 1;
    }

    // Read blocklist -- just like sdf_read_blocklist() but with the ability
    // to parse as we go
    sdf_read_summary(h);

    found = 0;
    range_start = 0;
    print_block = 0;
    if (metadata) {
        if (variable_ids || nrange) {
            DBG_FLUSH();
            print_block = 1;
        } else
            DBG_PRINT_FLUSH();
        found = 1;
    }

    // Construct the metadata blocklist using the contents of the buffer
    for (i = 0; i < h->nblocks; i++) {
        if (print_block) printf("\n  Block number: %i\n", i+1);
        sdf_read_block_info(h);

        if (print_block) {
            found = 0;

            for (n = range_start; n < nrange; n++) {
                if ((i+1) < range_list[n].start)
                    break;
                if ((i+1) <= range_list[n].end) {
                    found = 1;
                    break;
                }
                range_start++;
            }

            if (!found && variable_ids) {
                last_id = variable_ids;
                while (last_id) {
                    if (!memcmp(h->current_block->id, last_id->id,
                            strlen(last_id->id)+1)) {
                        found = 1;
                        break;
                    }
                    last_id = last_id->next;
                }
            }
        }

        if (found)
            DBG_PRINT_FLUSH();
        else
            DBG_FLUSH();
    }

    if (print_block) DBG_FLUSH();

    free(h->buffer);
    h->buffer = NULL;
    h->current_block = h->blocklist;

#ifdef SDF_DEBUG
    if (metadata) DBG_PRINT_FLUSH();
#endif

    // Done reading metadata

    if (!contents) return 0;

    found = 1;
    range_start = 0;
    if (variable_ids || nrange) print_block = 1;

    b = h->current_block = h->blocklist;

    for (i=0; i<h->nblocks; i++) {
        b = h->current_block;
        if (print_block) {
            found = 0;

            for (n = range_start; n < nrange; n++) {
                if ((i+1) < range_list[n].start)
                    break;
                if ((i+1) <= range_list[n].end) {
                    found = 1;
                    break;
                }
                range_start++;
            }

            if (!found && variable_ids) {
                last_id = variable_ids;
                while (last_id) {
                    if (!memcmp(h->current_block->id, last_id->id,
                            strlen(last_id->id)+1)) {
                        found = 1;
                        break;
                    }
                    last_id = last_id->next;
                }
            }

            if (found) {
                h->print = 1;
                sdf_read_data(h);
#ifdef SDF_DEBUG
                DBG_PRINT_FLUSH();
#endif
                h->print = 0;
            }
        } else
            sdf_read_data(h);

        h->current_block = b->next;

        if (found) DBG_PRINT_FLUSH();
    }
#ifdef SDF_DEBUG_ALL
#ifdef SDF_DEBUG
    if (debug) DBG_PRINT_FLUSH();
#endif
#endif
    sdf_close(h);
    printf("\n");

#ifdef PARALLEL
    MPI_Finalize();
#endif

    return 0;
}
