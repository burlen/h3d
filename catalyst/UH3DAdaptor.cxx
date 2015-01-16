/*=========================================================================

  Program:   ParaView
  Module:    UH3DAdaptor.cxx

  Copyright (c) Kitware, Inc.
  All rights reserved.
  See Copyright.txt or http://www.paraview.org/HTML/Copyright.html for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "FortranAdaptorAPI.h"
#include "vtkCellData.h"
#include "vtkCellType.h"
#include "vtkCPDataDescription.h"
#include "vtkCPInputDataDescription.h"
#include "vtkCPProcessor.h"
#include "vtkDoubleArray.h"
#include "vtkFieldData.h"
#include "vtkFloatArray.h"
#include "vtkImageData.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkUnstructuredGrid.h"
#include "vtkRectilinearGrid.h"
#include "vtkMultiBlockDataSet.h"

extern "C" void createvtkuniformgrid_(int * istrx,int * iendx,int *jstry,int *jendy,int *kstrz,int *kendz)
{
   if(!ParaViewCoProcessing::GetCoProcessorData())
     {
     vtkGenericWarningMacro("Unable to access CoProcessorData.");
     return;
     }
   vtkMultiBlockDataSet* grid = vtkMultiBlockDataSet::New();
   grid->SetNumberOfBlocks(1);

   vtkImageData* img = vtkImageData::New();
   img->Initialize();

   grid->SetBlock(0, img);
   img->Delete();

   img->SetSpacing(1.0,1.0, 1.0);
   img->SetExtent(*istrx,*iendx,*jstry,*jendy,*kstrz,*kendz);
   img->SetOrigin(0.0, 0.0, 0.0);

   ParaViewCoProcessing::GetCoProcessorData()->GetInputDescriptionByName("input")->SetGrid(grid);
   grid->Delete();
}

extern "C" void adduniformgridscalar_(char *fieldname,int *fieldname_len,double *data,int *size)
{
  vtkDoubleArray *arr = vtkDoubleArray::New();
  vtkStdString name(fieldname,*fieldname_len);
  arr->SetName(name);
  arr->SetNumberOfComponents(1);
  arr->SetArray(data,*size,1);
  vtkMultiBlockDataSet *grid = vtkMultiBlockDataSet::SafeDownCast(ParaViewCoProcessing::GetCoProcessorData()
                               ->GetInputDescriptionByName("input")->GetGrid());
  vtkDataSet *dataset = vtkDataSet::SafeDownCast(grid->GetBlock(0));
  dataset->GetPointData()->AddArray(arr);
  arr->Delete();
}