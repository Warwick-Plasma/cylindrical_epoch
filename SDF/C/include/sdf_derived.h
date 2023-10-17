/*
 * SDF (Self-Describing Format) Software Library
 * Copyright (c) 2012-2016, SDF Development Team
 *
 * Distributed under the terms of the BSD 3-clause License.
 * See the LICENSE file for details.
 */

/**
   @internal
   @file sdf_derived.h

   @brief Declarations for the SDF C-library.
   @details Routines for reading and writing SDF files.
   @author Dr Keith Bennett
   @date 15/02/2014
*/

#ifndef _SDF_DERIVED_H_
#define _SDF_DERIVED_H_
#include <sdf.h>

#ifdef __cplusplus
extern "C" {
#endif

int sdf_add_derived_blocks(sdf_file_t *h);
int sdf_add_derived_blocks_final(sdf_file_t *h);

#ifdef __cplusplus
}
#endif

#endif
