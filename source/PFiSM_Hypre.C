// Wrapper around Hypre functions to copy data to GPU when using
// Hypre with GPU support. If Hypre was built with GPU support,
// HYPRE_USING_CUDA is defined.
//
#include "_hypre_struct_mv.h"

#include <stdio.h>
#include <dlfcn.h>
#include <iostream>

/*--------------------------------------------------------------------------
 * HYPRE_StructVectorGetBoxValues
 *--------------------------------------------------------------------------*/

HYPRE_Int
PFiSM_HYPRE_StructVectorGetBoxValues( HYPRE_StructVector  vector,
                                HYPRE_Int          *ilower,
                                HYPRE_Int          *iupper,
                                HYPRE_Complex      *values )
{
   //std::cout<<"PFiSM version of HYPRE_StructVectorGetBoxValues()...\n";
   hypre_Index   new_ilower;
   hypre_Index   new_iupper;
   hypre_Box    *new_value_box;

   HYPRE_Int     d;
   HYPRE_Real* device_values;
   HYPRE_Int     volume;
#ifdef HYPRE_USING_CUDA
   cudaPointerAttributes attr;
   cudaError_t ret;
   ret = cudaPointerGetAttributes(&attr,values);
   if ( ret == cudaErrorInvalidValue )
   {
      cudaGetLastError(); /* clear out previous API error */

      volume=1;
      for (d = 0; d < hypre_StructVectorNDim(vector); d++)
      {
         volume = volume*(1+iupper[d]-ilower[d]);
      }
      //printf("HYPRE_StructVectorGetBoxValues, alloc mem device...\n");
      device_values   = hypre_CTAlloc(HYPRE_Real, volume,HYPRE_MEMORY_DEVICE);
      hypre_Memcpy(device_values,values,volume*sizeof(double),
                HYPRE_MEMORY_DEVICE,HYPRE_MEMORY_HOST);
   }
   else
   {
      device_values = values;
   }
#else
   device_values = values;
#endif

   hypre_SetIndex(new_ilower, 0);
   hypre_SetIndex(new_iupper, 0);
   for (d = 0; d < hypre_StructVectorNDim(vector); d++)
   {
      hypre_IndexD(new_ilower, d) = ilower[d];
      hypre_IndexD(new_iupper, d) = iupper[d];
   }
   new_value_box = hypre_BoxCreate(hypre_StructVectorNDim(vector));
   hypre_BoxSetExtents(new_value_box, new_ilower, new_iupper);

   hypre_StructVectorSetBoxValues(vector, new_value_box, new_value_box,
                                  device_values, -1, -1, 0);

   hypre_BoxDestroy(new_value_box);
#ifdef HYPRE_USING_CUDA
   if ( ret == cudaErrorInvalidValue )
   {
      hypre_Memcpy(values,device_values,volume*sizeof(double),
                   HYPRE_MEMORY_HOST,HYPRE_MEMORY_DEVICE);
      hypre_TFree(device_values,HYPRE_MEMORY_DEVICE);
   }
#endif
   return hypre_error_flag;
}

HYPRE_Int
HYPRE_StructVectorGetBoxValues( HYPRE_StructVector  vector,
                                HYPRE_Int          *ilower,
                                HYPRE_Int          *iupper,
                                HYPRE_Complex      *values )
{
   printf("Call HYPRE_StructVectorGetBoxValues...\n");
   static HYPRE_Int (*my_HYPRE_StructVectorGetBoxValues)
                        ( HYPRE_StructVector  vector,
                          HYPRE_Int          *ilower,
                          HYPRE_Int          *iupper,
                          HYPRE_Complex*      values );
   char* error;
   if(!my_HYPRE_StructVectorGetBoxValues)
   {
      *(void**)(&my_HYPRE_StructVectorGetBoxValues) =
         dlsym(RTLD_NEXT, "HYPRE_StructVectorGetBoxValues");
         if((error = dlerror()) != NULL ) {
            fputs(error, stderr);
         }
   }
   my_HYPRE_StructVectorGetBoxValues( vector, ilower, iupper, values);

   return PFiSM_HYPRE_StructVectorGetBoxValues( vector, ilower, iupper, values);
}

/*--------------------------------------------------------------------------
 * HYPRE_StructVectorSetBoxValues
 *--------------------------------------------------------------------------*/

HYPRE_Int
PFiSM_HYPRE_StructVectorSetBoxValues( HYPRE_StructVector  vector,
                                HYPRE_Int          *ilower,
                                HYPRE_Int          *iupper,
                                HYPRE_Complex      *values )
{
   //std::cout<<"PFiSM_HYPRE_StructVectorSetBoxValues\n";
   hypre_Index   new_ilower;
   hypre_Index   new_iupper;
   hypre_Box    *new_value_box;

   HYPRE_Int     d;
   HYPRE_Int     volume;
   HYPRE_Real   *device_values;
#ifdef HYPRE_USING_CUDA
   cudaPointerAttributes attr;
   cudaError_t ret;
   ret = cudaPointerGetAttributes(&attr,values);
   if ( ret == cudaErrorInvalidValue )
   {
      cudaGetLastError(); /* clear out previous API error */

      volume=1;
      for (d = 0; d < hypre_StructVectorNDim(vector); d++)
      {
         volume = volume*(1+iupper[d]-ilower[d]);
      }
      //printf("HYPRE_StructVectorSetBoxValues, alloc mem device...\n");
       device_values   = hypre_CTAlloc(HYPRE_Real, volume,HYPRE_MEMORY_DEVICE);
      hypre_Memcpy(device_values,values,volume*sizeof(double),
                HYPRE_MEMORY_DEVICE,HYPRE_MEMORY_HOST);
   }
   else
   {
      device_values = values;
   }
#else
   device_values = values;
#endif

   hypre_SetIndex(new_ilower, 0);
   hypre_SetIndex(new_iupper, 0);
   for (d = 0; d < hypre_StructVectorNDim(vector); d++)
   {
      hypre_IndexD(new_ilower, d) = ilower[d];
      hypre_IndexD(new_iupper, d) = iupper[d];
   }
   new_value_box = hypre_BoxCreate(hypre_StructVectorNDim(vector));
   hypre_BoxSetExtents(new_value_box, new_ilower, new_iupper);

   hypre_StructVectorSetBoxValues(vector, new_value_box, new_value_box,
                                  device_values, 0, -1, 0);

   hypre_BoxDestroy(new_value_box);

#ifdef HYPRE_USING_CUDA
   if ( ret == cudaErrorInvalidValue )
   {
       hypre_TFree(device_values,HYPRE_MEMORY_DEVICE);
   }
#endif

   return hypre_error_flag;
}

HYPRE_Int
HYPRE_StructVectorSetBoxValues( HYPRE_StructVector  vector,
                                HYPRE_Int          *ilower,
                                HYPRE_Int          *iupper,
                                HYPRE_Complex      *values )
{
   printf("Call HYPRE_StructVectorSetBoxValues...\n");
   static HYPRE_Int (*my_HYPRE_StructVectorSetBoxValues)
                        ( HYPRE_StructVector  vector,
                          HYPRE_Int          *ilower,
                          HYPRE_Int          *iupper,
                          HYPRE_Complex*      values );
   char* error;
   if(!my_HYPRE_StructVectorSetBoxValues)
   {
      *(void**)(&my_HYPRE_StructVectorSetBoxValues) =
         dlsym(RTLD_NEXT, "HYPRE_StructVectorSetBoxValues");
         if((error = dlerror()) != NULL ) {
            fputs(error, stderr);
         }
   }
   my_HYPRE_StructVectorSetBoxValues( vector, ilower, iupper, values);

   return PFiSM_HYPRE_StructVectorSetBoxValues( vector, ilower, iupper, values);
}

/*--------------------------------------------------------------------------
 * HYPRE_StructMatrixSetBoxValues
 *--------------------------------------------------------------------------*/

HYPRE_Int
PFiSM_HYPRE_StructMatrixSetBoxValues( HYPRE_StructMatrix  matrix,
                                HYPRE_Int          *ilower,
                                HYPRE_Int          *iupper,
                                HYPRE_Int           num_stencil_indices,
                                HYPRE_Int          *stencil_indices,
                                HYPRE_Complex      *values )
{
   hypre_Index         new_ilower;
   hypre_Index         new_iupper;
   hypre_Box          *new_value_box;
   HYPRE_Int           d;
   HYPRE_Int     volume;
   HYPRE_Real   *device_values;
#ifdef HYPRE_USING_CUDA
   cudaPointerAttributes attr;
   cudaError_t ret;
   ret = cudaPointerGetAttributes(&attr,values);
   if ( ret == cudaErrorInvalidValue )
   {
      cudaGetLastError(); /* clear out previous API error */
      volume=1;
      for (d = 0; d < hypre_StructMatrixNDim(matrix); d++)
      {
         volume = volume*(1+iupper[d]-ilower[d]);
      }
      volume = volume*num_stencil_indices;
      //printf("HYPRE_StructMatrixSetBoxValues, alloc mem device...\n");
      device_values   = hypre_CTAlloc(HYPRE_Real, volume,HYPRE_MEMORY_DEVICE);
      hypre_Memcpy(device_values,values,volume*sizeof(double),
                HYPRE_MEMORY_DEVICE,HYPRE_MEMORY_HOST);
   }
   else
   {
      device_values = values;
   }
#else
   device_values = values;
#endif

   hypre_SetIndex(new_ilower, 0);
   hypre_SetIndex(new_iupper, 0);
   for (d = 0; d < hypre_StructMatrixNDim(matrix); d++)
   {
      hypre_IndexD(new_ilower, d) = ilower[d];
      hypre_IndexD(new_iupper, d) = iupper[d];
   }
   new_value_box = hypre_BoxCreate(hypre_StructMatrixNDim(matrix));
   hypre_BoxSetExtents(new_value_box, new_ilower, new_iupper);

   hypre_StructMatrixSetBoxValues(matrix, new_value_box, new_value_box,
                                  num_stencil_indices, stencil_indices,
                                  device_values, 0, -1, 0);

   hypre_BoxDestroy(new_value_box);
#ifdef HYPRE_USING_CUDA
   if ( ret == cudaErrorInvalidValue )
   {
       hypre_TFree(device_values,HYPRE_MEMORY_DEVICE);
   }
#endif

   return hypre_error_flag;
}
