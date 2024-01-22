/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkComputeRestrictedHistogram.txx,v $
  Language:  C++
  Date:      $Date: 2006/01/11 19:43:31 $
  Version:   $Revision: 1.14 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itk_MultiStencilEuclideanFastMarching1stOrderImageFilter_txx
#define _itk_MultiStencilEuclideanFastMarching1stOrderImageFilter_txx
#include "itkMultiStencilEuclideanFastMarching1stOrderImageFilter.h"

#include "itkNumericTraits.h"
#include "vnl/vnl_math.h"
#include "itkMath.h"

#include "itkConstNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"
#include "itkNeighborhoodAlgorithm.h"

#include <iostream>

///PARA DEPURAR
//#include "itkImageFileWriter.h"


namespace itk
{


template< class TInputImage, class TRealImage, class TOutputImage>
MultiStencilEuclideanFastMarching1stOrderImageFilter< TInputImage, TRealImage, TOutputImage>
::MultiStencilEuclideanFastMarching1stOrderImageFilter()
{
}

template< class TInputImage, class TRealImage, class TOutputImage>
bool 
MultiStencilEuclideanFastMarching1stOrderImageFilter< TInputImage, TRealImage, TOutputImage>
::ComputeNodeValue( const InputIndexType& idx )
{
  //FIXME: ELIMINAR CUANDO ESTÉ DEPURADO
  InputIndexType idxToDebug;
  idxToDebug[0] = 2; idxToDebug[1] = 2; idxToDebug[2] = 3; 
  if (idx==idxToDebug)
    std::cout << "Calculando IDX de interés" << std::endl;
  
  //Image pointers
  typename OutputImageType::Pointer T = this->GetOutput();


  //Iterate on the stencils
  DoubleContainerType solTvec; solTvec.clear();
  for (unsigned int stencil=0; stencil<Superclass::m_NumberOfStencils; stencil++)
  {
    UpwindVectorType upwind;
    SignVectorType sign;

    bool isValid = this->ComputeStencilUpwindNodes(stencil,idx,upwind,sign);
    if (!isValid)
      continue;
    
    bool notSolved = true;
    do
    {
      double cost = static_cast<double>(Superclass::m_CostImage->GetPixel(idx));
      if (cost < Superclass::m_NumericTol)
	cost = Superclass::m_NumericTol;
      
      //Get coefficients for the linear equation 
      double f0, f1, f2;
      this->Superclass::ComputeStencilEquationCoefficients( stencil, 
							    upwind,sign,Superclass::m_StencilOffsetNorms[stencil],cost,
							    &f0,&f1,&f2 );
      
      //Solve the equation
      double candT;
      bool isReal = this->Superclass::SolveQuadraticEquation(f0,f1,f2, &candT);
      
      if (isReal)
      {
	//Check the solution and assign it if valid. If not, try with less nodes
	std::vector<unsigned int> lowerTidx; lowerTidx.clear();
	
	for (unsigned int i=0;i<InputDimension;i++)
	{
	  if ((sign[i]!=0)&&(candT < upwind[i] - Superclass::m_NumericTol))
	    lowerTidx.push_back(i);
	}
	
	if (lowerTidx.empty())
	{
	  //The solution is a valid candidate
	  solTvec.push_back( candT );
	  notSolved = false;
	}
	else if ( InputDimension-lowerTidx.size() > 0 ) 
	{
	  for (unsigned int i=0; i<lowerTidx.size(); i++)
	  {
	    upwind[lowerTidx[i]] = NumericTraits<OutputPixelType>::NonpositiveMin();
	    sign[lowerTidx[i]] = 0;
	  }
	  
	  //notSolved = true; //But it is not necessary to state it
	}
	else
	{
	  //Cannot solve it. Continuing to the next stencil
	  notSolved = false;
	}
      }
      else
      {
	//TCAND(s) = inf;
	notSolved = false;
      }

    } while(notSolved);
    
  }
  
  //Select the stencil with the minimum Tm 
  //and assign the solution to the image outputs
  if ( solTvec.empty() )
  {
    //There are no adequate solutions for this node at the moment
    return false;
  }
  else
  {
    //Select the Tm and Ts solutions from the 
    //stencil with minimum Tm. If two Tm are equal,
    //we choose the one with least Ts.
    double minValue = NumericTraits<double>::max();
    unsigned int minPos=0;
    for (unsigned int i=0; i<solTvec.size(); i++)
    {
      if (minValue>solTvec[i])
      {
	minValue = solTvec[i];
	minPos = i;
      }
      else if (minValue==solTvec[i])
      {
	minPos = i;
      }
    }
    
    //Assign them to the outputs
    OutputPixelType chosenT = (OutputPixelType)solTvec.at(minPos);
    T->SetPixel( idx, chosenT );
    this->Superclass::m_ChosenStencilImage->SetPixel( idx, minPos );
    
    return true;
  }
}




template< class TInputImage, class TRealImage, class TOutputImage>
bool 
MultiStencilEuclideanFastMarching1stOrderImageFilter< TInputImage, TRealImage, TOutputImage>
::ComputeStencilUpwindNodes( unsigned int stcl,
			     const InputIndexType& idx,
			     UpwindVectorType &upwind,
			     SignVectorType &sign )
{
  // WARNING: IMPLEMENTADA, SIN DEPURAR
  // Const pointers to images
  typename OutputImageType::ConstPointer output = this->GetOutput();
  
  // Initialize upwind and sign vectors
  sign.Fill(0);
  upwind.Fill(NumericTraits<OutputPixelType>::NonpositiveMin());
  
  // Iterate through the stencil dimension
  OutputPixelType candval[2]; 
  InputIndexType fwIdx,bwIdx;
  bool hasNeighs = false;
  
  for (unsigned int dim=0; dim<InputDimension; dim++)
  {
    //Initialize upwind candidates
    for (unsigned int i=0; i<2; i++)
      candval[i] = Superclass::m_LargeValue;
    
    //Check if the upwind node is within the image and frozen
    fwIdx = idx + Superclass::m_StencilOffsets[stcl][dim];
    if ((Superclass::m_LargestRegion.IsInside(fwIdx)) && (Superclass::m_LabelImage->GetPixel(fwIdx)==Superclass::KnownPoint))
    {
      candval[0] = output->GetPixel(fwIdx);
    }
    
    //Check if the downwind node is within the image and frozen
    bwIdx = idx - Superclass::m_StencilOffsets[stcl][dim];
    if ((Superclass::m_LargestRegion.IsInside(bwIdx)) && (Superclass::m_LabelImage->GetPixel(bwIdx)==Superclass::KnownPoint))
    {
      candval[1] = output->GetPixel(bwIdx);
    }
    
    //Depending on candval, assign upwind[dim] and sign[dim]
    if ( (candval[0]==Superclass::m_LargeValue) && (candval[1]==Superclass::m_LargeValue) )
    {
      continue;
    }
    else if (candval[0]<candval[1])
    {
      sign[dim] = 1;
      upwind[dim] = candval[0];
      hasNeighs = true;
    }
    else
    {
      sign[dim] = -1;
      upwind[dim] = candval[1];
      hasNeighs = true;
    }
  }
  
  return hasNeighs;
}



} // end namespace itk

#endif //#ifndef





