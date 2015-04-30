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

#include "vtkCellData.h"
#include "vtkCellType.h"
#include "vtkCPDataDescription.h"
#include "vtkCPExodusIINodalCoordinatesTemplate.h"
#include "vtkCPInputDataDescription.h"
#include "vtkCPProcessor.h"
#include "vtkCPPythonAdaptorAPI.h"
#include "vtkDoubleArray.h"
#include "vtkFieldData.h"
#include "vtkFloatArray.h"
#include "vtkImageData.h"
#include "vtkPointData.h"

extern "C"
void createvtkuniformgrid_(
  int * istrx,int * iendx,int *jstry,int *jendy,int *kstrz,int *kendz)
{
  if(!vtkCPPythonAdaptorAPI::GetCoProcessorData())
    {
    vtkGenericWarningMacro("Unable to access CoProcessorData.");
    return;
    }

  vtkImageData* img = vtkImageData::New();
  img->Initialize();

  img->SetSpacing(1.0,1.0, 1.0);
  img->SetExtent(*istrx,*iendx,*jstry,*jendy,*kstrz,*kendz);
  img->SetOrigin(0.0, 0.0, 0.0);

  vtkCPPythonAdaptorAPI::GetCoProcessorData()->GetInputDescriptionByName("input")->SetGrid(img);
  img->Delete();
}

extern "C"
void adduniformgridscalar_(
  const char *fieldname,int *fieldname_len,double *data,int *size)
{
  vtkDoubleArray *arr = vtkDoubleArray::New();
  vtkStdString name(fieldname,*fieldname_len);
  arr->SetName(name);
  arr->SetNumberOfComponents(1);
  arr->SetArray(data,*size,1);
  vtkDataSet *dataset = vtkDataSet::SafeDownCast(
    vtkCPPythonAdaptorAPI::GetCoProcessorData()->GetInputDescriptionByName("input")->GetGrid());
  dataset->GetPointData()->AddArray(arr);
  arr->Delete();
}

extern "C"
void adduniformgridvector_(
  const char *fieldname,int *fieldname_len,
  double *datax, double *datay, double *dataz, int *size)
{
  vtkStdString name(fieldname, *fieldname_len);
  int newLength = *fieldname_len + 1;
  //Still add the scalars for debugging purposes
  adduniformgridscalar_( (name + "x").c_str(), &newLength, datax, size);
  adduniformgridscalar_( (name + "y").c_str(), &newLength, datay, size);
  adduniformgridscalar_( (name + "z").c_str(), &newLength, dataz, size);

/*
  // Add the 'zero copy' vector arrays
  vtkCPExodusIINodalCoordinatesTemplate<double> *a = 
    vtkCPExodusIINodalCoordinatesTemplate<double>::New();
  a->SetName((name + "0").c_str());
  a->SetNumberOfComponents(3);
  a->SetExodusScalarArrays(datax, datay, dataz, *size);
  vtkMultiBlockDataSet *grid = vtkMultiBlockDataSet::SafeDownCast(
    vtkCPPythonAdaptorAPI::GetCoProcessorData()->GetInputDescriptionByName("input")->GetGrid());
  vtkDataSet *dataset = vtkDataSet::SafeDownCast(grid->GetBlock(0));
  dataset->GetPointData()->AddArray(a);
  a->Delete();
*/
}
