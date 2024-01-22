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
#ifndef _itk_MultiStencilEuclideanFastMarchingCenteredImageFilter_txx
#define _itk_MultiStencilEuclideanFastMarchingCenteredImageFilter_txx
#include "itkMultiStencilEuclideanFastMarchingCenteredImageFilter.h"

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
MultiStencilEuclideanFastMarchingCenteredImageFilter< TInputImage, TRealImage, TOutputImage>
::MultiStencilEuclideanFastMarchingCenteredImageFilter()
{
}

template< class TInputImage, class TRealImage, class TOutputImage>
bool 
MultiStencilEuclideanFastMarchingCenteredImageFilter< TInputImage, TRealImage, TOutputImage>
::ComputeNodeValue( const InputIndexType& idx )
{
  //FIXME: ELIMINAR CUANDO ESTÉ DEPURADO
  InputIndexType idxToDebug;
  idxToDebug[0] = 3; idxToDebug[1] = 24; idxToDebug[2] = 1; 
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
      UpwindVectorType eqUpwind;
      double eqCost;
      DoubleContainerType eqNorms;
      
      if (InputDimension==2)
	this->ComputeEquivalentValues2D(stencil,idx,upwind,sign,eqUpwind,&eqCost,eqNorms);
      else      
	this->ComputeEquivalentValues3D(stencil,idx,upwind,sign,eqUpwind,&eqCost,eqNorms);
            
      
      //Get coefficients for the linear equation %TODO
      double f0, f1, f2;
      this->Superclass::ComputeStencilEquationCoefficients( stencil, eqUpwind,sign,eqNorms,eqCost, &f0,&f1,&f2 );
      
      //Solve the equation
      double candT;
      bool isReal = this->Superclass::SolveQuadraticEquation(f0,f1,f2, &candT);
      
      if (isReal)
      {
	//Check the solution and assign it if valid. If not, try with less nodes
	std::vector<unsigned int> lowerTidx; lowerTidx.clear();
	
	for (unsigned int i=0;i<InputDimension;i++)
	{
	  if ((sign[i]!=0)&&(candT < eqUpwind[i] - Superclass::m_NumericTol))
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
MultiStencilEuclideanFastMarchingCenteredImageFilter< TInputImage, TRealImage, TOutputImage>
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


template< class TInputImage, class TRealImage, class TOutputImage>
void 
MultiStencilEuclideanFastMarchingCenteredImageFilter< TInputImage, TRealImage, TOutputImage>
::ComputeEquivalentValues3D( unsigned int stcl,
				const InputIndexType& idx,
				const UpwindVectorType &upwind,
				const SignVectorType& sign,
				UpwindVectorType &eqUpwind,
				double *eqCost,
				DoubleContainerType &eqNorms )
{
  // Const pointers to images
  typename OutputImageType::ConstPointer output = this->GetOutput();
//   typename RealImageType::ConstPointer cost = this->GetCostInput();
  
  // Apply the sign to the offsets and compute the minimum valid upwind value
  OffsetContainerType validOffsets; 
  std::vector<unsigned int> validIdx;
  OutputPixelType minUpwind = Superclass::m_LargeValue;
  
  validOffsets.clear(); validOffsets.reserve( 16 ); //To avoid resizing the vector
  validIdx.clear(); validIdx.reserve( 16 );
  for (unsigned int i=0; i<InputDimension;i++)
  {
    if (sign[i]!=0)
    {
      OffsetType signedOff;
      for (unsigned j=0;j<InputDimension;j++)
	signedOff[j] = sign[i] * Superclass::m_StencilOffsets[stcl][i][j];
      
      validOffsets.push_back( signedOff );
      validIdx.push_back(i);
      
      if (minUpwind>upwind[i])
	minUpwind = upwind[i];
    }
  }
  
  // Initialize outputs
  eqNorms = Superclass::m_StencilOffsetNorms[stcl];
  eqUpwind = upwind;
  
  //Depending on how many frozen neighbors there are, the processing changes.
  unsigned int numNeighs = validIdx.size();
  typename Superclass::ContinuousIndexType cIdx;
  
  if (numNeighs==1)
  {
    for (unsigned int i=0; i<InputDimension;i++)
      cIdx[i] = idx[i] + static_cast<double>( validOffsets[0][i] )/2;
    
    *eqCost = (double) Superclass::m_Interpolator->EvaluateAtContinuousIndex( cIdx );
  }
  else if (numNeighs==2)
  {
    InputIndexType diagIdx = idx + validOffsets[0] + validOffsets[1];
    
    if (Superclass::m_LargestRegion.IsInside(diagIdx) && (output->GetPixel(diagIdx)<minUpwind))
    {
      for (unsigned int i=0; i<InputDimension;i++)
	cIdx[i] = static_cast<double>(diagIdx[i] + idx[i])/2;
      
      *eqCost = (double) Superclass::m_Interpolator->EvaluateAtContinuousIndex( cIdx );
      
      for (unsigned int i=0; i<InputDimension;i++)
	eqNorms[i] *= 2;
      
      //Equivalent values
      OutputPixelType Tdiag = output->GetPixel(diagIdx);
      eqUpwind[ validIdx[0] ] = Tdiag + upwind[validIdx[0]] -upwind[validIdx[1]];
      eqUpwind[ validIdx[1] ] = Tdiag + upwind[validIdx[1]] -upwind[validIdx[0]];
    }
    else
    {
      //Revert to first order scheme
      *eqCost = static_cast<double>( Superclass::m_CostImage->GetPixel(idx) );
    }
  }
  else
  {
    InputIndexType diagIdx = idx + validOffsets[0] + validOffsets[1] + validOffsets[2];
    IndexContainerType diagFaceIdxs; 
    diagFaceIdxs.clear(); diagFaceIdxs.reserve( 4 );
    diagFaceIdxs.push_back( idx + validOffsets[0] + validOffsets[1] );
    diagFaceIdxs.push_back( idx + validOffsets[0] + validOffsets[2] );
    diagFaceIdxs.push_back( idx + validOffsets[1] + validOffsets[2] );
    
    bool allFrozen = true;
    for (unsigned int i=0; i<diagFaceIdxs.size();i++)
    {
      if ( (!Superclass::m_LargestRegion.IsInside(diagFaceIdxs[i])) || (Superclass::m_LabelImage->GetPixel(diagFaceIdxs[i])!=Superclass::KnownPoint))
      {
	allFrozen = false;
	break;
      }
    }
    
    if ( allFrozen && Superclass::m_LargestRegion.IsInside(diagIdx) && (output->GetPixel(diagIdx)<minUpwind))
    {
      for (unsigned int i=0; i<InputDimension;i++)
	cIdx[i] = static_cast<double>(diagIdx[i] + idx[i])/2;
      
      *eqCost = (double) Superclass::m_Interpolator->EvaluateAtContinuousIndex( cIdx );
      
      for (unsigned int i=0; i<InputDimension;i++)
	eqNorms[i] *= 4;
      
      //values in the diagonals
      OutputPixelType Tdiag = output->GetPixel(diagIdx);
      
      OutputPixelType Tfaces[3];
      for (unsigned int i=0; i<3; i++)
	Tfaces[i] = output->GetPixel(diagFaceIdxs[i]);
      
      eqUpwind[0] = Tdiag + Tfaces[0] + Tfaces[1] + upwind[0] - upwind[1] - upwind[2] - Tfaces[2];
      eqUpwind[1] = Tdiag + Tfaces[0] + Tfaces[2] + upwind[1] - upwind[0] - upwind[2] - Tfaces[1];
      eqUpwind[2] = Tdiag + Tfaces[1] + Tfaces[2] + upwind[2] - upwind[0] - upwind[1] - Tfaces[0];
    }
    else
    {
      //Revert to first order scheme
      *eqCost = Superclass::m_CostImage->GetPixel(idx);
    }
  }
  
  // The FM method is defined for positive costs. Due to this,
  // if the cost is null we set it to the numeric tolerance.
  if (*eqCost < Superclass::m_NumericTol)
    *eqCost = Superclass::m_NumericTol;
}



template< class TInputImage, class TRealImage, class TOutputImage>
void 
MultiStencilEuclideanFastMarchingCenteredImageFilter< TInputImage, TRealImage, TOutputImage>
::ComputeEquivalentValues2D( unsigned int stcl,
				const InputIndexType& idx,
				const UpwindVectorType &upwind,
				const SignVectorType& sign,
				UpwindVectorType &eqUpwind,
				double *eqCost,
				DoubleContainerType &eqNorms )
{
  // Const pointers to images
  typename OutputImageType::ConstPointer output = this->GetOutput();
//   typename RealImageType::ConstPointer cost = this->GetCostInput();
  
  // Apply the sign to the offsets and compute the minimum valid upwind value
  OffsetContainerType validOffsets; 
  std::vector<unsigned int> validIdx;
  OutputPixelType minUpwind = Superclass::m_LargeValue;
  
  validOffsets.clear(); validOffsets.reserve( 16 ); //To avoid resizing the vector
  validIdx.clear(); validIdx.reserve( 16 );
  for (unsigned int i=0; i<InputDimension;i++)
  {
    if (sign[i]!=0)
    {
      OffsetType signedOff;
      for (unsigned j=0;j<InputDimension;j++)
	signedOff[j] = sign[i] * Superclass::m_StencilOffsets[stcl][i][j];
      
      validOffsets.push_back( signedOff );
      validIdx.push_back(i);
      
      if (minUpwind>upwind[i])
	minUpwind = upwind[i];
    }
  }
  
  // Initialize outputs
  eqNorms = Superclass::m_StencilOffsetNorms[stcl];
  eqUpwind = upwind;
  
  //Depending on how many frozen neighbors there are, the processing changes.
  unsigned int numNeighs = validIdx.size();
  typename Superclass::ContinuousIndexType cIdx;
  
  if (numNeighs==1)
  {
    for (unsigned int i=0; i<InputDimension;i++)
      cIdx[i] = idx[i] + static_cast<double>( validOffsets[0][i] )/2;
    
    *eqCost = (double) Superclass::m_Interpolator->EvaluateAtContinuousIndex( cIdx );
  }
  else
  {
    InputIndexType diagIdx = idx + validOffsets[0] + validOffsets[1];
    
    if (Superclass::m_LargestRegion.IsInside(diagIdx) && (output->GetPixel(diagIdx)<minUpwind))
    {
      for (unsigned int i=0; i<InputDimension;i++)
	cIdx[i] = static_cast<double>(diagIdx[i] + idx[i])/2;
      
      *eqCost = (double) Superclass::m_Interpolator->EvaluateAtContinuousIndex( cIdx );
      
      for (unsigned int i=0; i<InputDimension;i++)
	eqNorms[i] *= 2;
      
      //Equivalent values
      OutputPixelType Tdiag = output->GetPixel(diagIdx);
      eqUpwind[0] = Tdiag + upwind[0] -upwind[1];
      eqUpwind[1] = Tdiag + upwind[1] -upwind[0];
    }
    else
    {
      //Revert to first order scheme
      *eqCost = static_cast<double>( Superclass::m_CostImage->GetPixel(idx) );
    }
  }
  
  // The FM method is defined for positive costs. Due to this,
  // if the cost is null we set it to the numeric tolerance.
  if (*eqCost < Superclass::m_NumericTol)
    *eqCost = Superclass::m_NumericTol;
}



} // end namespace itk

#endif //#ifndef





