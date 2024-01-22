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
#ifndef _itk_MultiStencilEuclideanFastMarching2ndOrderImageFilter_txx
#define _itk_MultiStencilEuclideanFastMarching2ndOrderImageFilter_txx
#include "itkMultiStencilEuclideanFastMarching2ndOrderImageFilter.h"

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
MultiStencilEuclideanFastMarching2ndOrderImageFilter< TInputImage, TRealImage, TOutputImage>
::MultiStencilEuclideanFastMarching2ndOrderImageFilter()
{
}

template< class TInputImage, class TRealImage, class TOutputImage>
bool 
MultiStencilEuclideanFastMarching2ndOrderImageFilter< TInputImage, TRealImage, TOutputImage>
::ComputeNodeValue( const InputIndexType& idx )
{
  //FIXME: ELIMINAR CUANDO ESTÉ DEPURADO
  InputIndexType idxToDebug;
  idxToDebug[0] = 3; idxToDebug[1] = 1; idxToDebug[2] = 3; 
  if (idx==idxToDebug)
    std::cout << "Calculando IDX de interés" << std::endl;
  
  //Image pointers
  typename OutputImageType::Pointer T = this->GetOutput();


  //Iterate on the stencils
  DoubleContainerType solTvec; solTvec.clear();
  for (unsigned int stencil=0; stencil<Superclass::m_NumberOfStencils; stencil++)
  {
    UpwindVectorType upwind1, upwind2;
    SignVectorType sign, secord;

    bool isValid = this->ComputeStencilUpwindNodes(stencil,idx,upwind1,upwind2,sign,secord);
    if (!isValid)
      continue;
    
    bool notSolved = true;
    do
    {
      UpwindVectorType eqUpwind;
      double eqCost;
      DoubleContainerType eqNorms;
      
      this->ComputeEquivalentValues(stencil,idx,upwind1,upwind2,sign,secord,eqUpwind,&eqCost,eqNorms);
      
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
	    upwind1[lowerTidx[i]] = NumericTraits<OutputPixelType>::NonpositiveMin();
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
MultiStencilEuclideanFastMarching2ndOrderImageFilter< TInputImage, TRealImage, TOutputImage>
::ComputeStencilUpwindNodes( unsigned int stcl,
			     const InputIndexType& idx,
			     UpwindVectorType &upwind1,
			     UpwindVectorType &upwind2,
			     SignVectorType &sign,
			     SignVectorType &secord )
{
  // WARNING: IMPLEMENTADA, SIN DEPURAR
  // Const pointers to images
  typename OutputImageType::ConstPointer output = this->GetOutput();
  
  // Initialize upwind and sign vectors
  sign.Fill(0);
  secord.Fill(0);
  upwind1.Fill(NumericTraits<OutputPixelType>::NonpositiveMin());
  upwind2.Fill(NumericTraits<OutputPixelType>::NonpositiveMin());
  
  // Iterate through the stencil dimension
  OutputPixelType candval1[2], candval2[2]; 
  InputIndexType fwIdx,bwIdx, fw2Idx,bw2Idx;
  bool hasNeighs = false;
  
  for (unsigned int dim=0; dim<InputDimension; dim++)
  {
    //Initialize upwind candidates
    for (unsigned int i=0; i<2; i++)
    {
      candval1[i] = Superclass::m_LargeValue;
      candval2[i] = Superclass::m_LargeValue;
    }
    
    //Check if the upwind node is within the image and frozen
    fwIdx = idx + Superclass::m_StencilOffsets[stcl][dim];
    if ((Superclass::m_LargestRegion.IsInside(fwIdx)) && (Superclass::m_LabelImage->GetPixel(fwIdx)==Superclass::KnownPoint))
    {
      candval1[0] = output->GetPixel(fwIdx);
      
      //Check the node at two forward steps
      fw2Idx = fwIdx + Superclass::m_StencilOffsets[stcl][dim];
      if ((Superclass::m_LargestRegion.IsInside(fw2Idx)) && (Superclass::m_LabelImage->GetPixel(fw2Idx)==Superclass::KnownPoint))
      {
	candval2[0] = output->GetPixel(fw2Idx);
      }
    }
    
    //Check if the downwind node is within the image and frozen
    bwIdx = idx - Superclass::m_StencilOffsets[stcl][dim];
    if ((Superclass::m_LargestRegion.IsInside(bwIdx)) && (Superclass::m_LabelImage->GetPixel(bwIdx)==Superclass::KnownPoint))
    {
      candval1[1] = output->GetPixel(bwIdx);
      
      //Check the node at two backward steps
      bw2Idx = bwIdx - Superclass::m_StencilOffsets[stcl][dim];
      if ((Superclass::m_LargestRegion.IsInside(bw2Idx)) && (Superclass::m_LabelImage->GetPixel(bw2Idx)==Superclass::KnownPoint))
      {
	candval2[1] = output->GetPixel(bw2Idx);
      }
    }
    
    //Depending on candval, assign upwind[dim] and sign[dim]
    if ( (candval1[0]==Superclass::m_LargeValue) && (candval1[1]==Superclass::m_LargeValue) )
    {
      continue;
    }
    else if (candval1[0]<candval1[1])
    {
      sign[dim] = 1;
      upwind1[dim] = candval1[0];
      hasNeighs = true;
      if ((candval2[0]!=Superclass::m_LargeValue)&&(candval2[0]<upwind1[dim]))
      {
	upwind2[dim] = candval2[0];
	secord[dim] = 1;
      }
    }
    else
    {
      sign[dim] = -1;
      upwind1[dim] = candval1[1];
      hasNeighs = true;
      if ((candval2[1]!=Superclass::m_LargeValue)&&(candval2[1]<upwind1[dim]))
      {
	upwind2[dim] = candval2[1];
	secord[dim] = 1;
      }
    }
  }
  
  return hasNeighs;
}


template< class TInputImage, class TRealImage, class TOutputImage>
void 
MultiStencilEuclideanFastMarching2ndOrderImageFilter< TInputImage, TRealImage, TOutputImage>
::ComputeEquivalentValues( unsigned int stcl,
				const InputIndexType& idx,
				const UpwindVectorType &upwind1,
				const UpwindVectorType &upwind2,
				const SignVectorType& sign,
				const SignVectorType& secord,
				UpwindVectorType &eqUpwind,
				double *eqCost,
				DoubleContainerType &eqNorms )
{
  // Const pointers to images
  typename OutputImageType::ConstPointer output = this->GetOutput();
  
  
  // Initialize or assign outputs
  eqNorms = Superclass::m_StencilOffsetNorms[stcl];
  eqUpwind = upwind1;
  *eqCost = static_cast<double>( Superclass::m_CostImage->GetPixel(idx) );
  
  // Iterate over the dimensions. If there is a valid second upwind node,
  // compute the equivalent values for that stencil vector.
  for (unsigned int i=0;i<InputDimension;i++)
  {
    if (secord[i]>0)
    {
      eqUpwind[i] = (4*upwind1[i]-upwind2[i]) /3;
      eqNorms[i] = 2*eqNorms[i]/3.0f;
    }
  }
  
  // The FM method is defined for positive costs. Due to this,
  // if the cost is null we set it to the numeric tolerance.
  if (*eqCost < Superclass::m_NumericTol)
    *eqCost = Superclass::m_NumericTol;
}





} // end namespace itk

#endif //#ifndef





