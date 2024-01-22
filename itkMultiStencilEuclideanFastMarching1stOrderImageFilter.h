/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkComputeRestrictedHistogram.h,v $
  Language:  C++
  Date:      $Date: 2006/03/29 14:53:40 $
  Version:   $Revision: 1.5 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkMultiStencilEuclideanFastMarching1stOrderImageFilter_h
#define _itkMultiStencilEuclideanFastMarching1stOrderImageFilter_h

#include "itkMultiStencilEuclideanFastMarchingImageFilterBase.h"
#include "itkImage.h"
#include "itkArray2D.h"
#include "itkArray.h"
#include "itkIndex.h"
#include "itkGradientImageFilter.h"
#include "itkMatrix.h"
#include "itkFixedArray.h"
#include "itkImageRegion.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkMacro.h"

#include "itkBinaryHeap.h"

#include "itkLevelSetNode.h"
#include <functional>
#include <queue>
#include <ctime>
// #include <cstring>
#include <string.h>


namespace itk
{
/** \class MultiStencilEuclideanFastMarching1stOrderImageFilter
 * \brief Implementation of the configurable multistencil Fast Marching method.
 * To write
 */
template <class TInputImage, class TRealImage, class TOutputImage=TRealImage>
class ITK_EXPORT MultiStencilEuclideanFastMarching1stOrderImageFilter
		  : public MultiStencilEuclideanFastMarchingImageFilterBase< TInputImage, TRealImage, TOutputImage >
{
public:
  /** Convenient typedefs for simplifying declarations. */
  typedef TInputImage                                          InputImageType;
  typedef TRealImage                                           RealImageType;
  typedef TOutputImage                                         OutputImageType;

  /** Standard class typedefs. */
  typedef MultiStencilEuclideanFastMarching1stOrderImageFilter 		Self;
  typedef MultiStencilEuclideanFastMarchingImageFilterBase< 
		InputImageType,RealImageType, OutputImageType> 		Superclass;
  typedef SmartPointer<Self>                                   		Pointer;
  typedef SmartPointer<const Self>                             		ConstPointer;

  /** Template image dimensions */
  static const unsigned int InputDimension	= TInputImage::ImageDimension;
  static const unsigned int RealDimension	= TRealImage::ImageDimension;
  static const unsigned int OutputDimension	= TOutputImage::ImageDimension;
  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MultiStencilEuclideanFastMarching1stOrderImageFilter, MultiStencilEuclideanFastMarchingImageFilterBase );

  /** Image typedef support. */
  typedef typename InputImageType::PixelType               InputPixelType;
  typedef typename RealImageType::PixelType                RealPixelType;
  typedef typename OutputImageType::PixelType              OutputPixelType;
  typedef typename NumericTraits<InputPixelType>::RealType InputRealType;
  typedef typename InputImageType::RegionType              InputImageRegionType;
  typedef typename OutputImageType::RegionType             OutputImageRegionType;
  typedef typename InputImageType::SizeType                InputSizeType;
  typedef typename InputImageType::SpacingType		   InputSpacingType;
  typedef typename RealImageType::SizeType                 RealSizeType;
  typedef typename InputImageType::IndexType               InputIndexType;
  typedef typename InputImageType::ConstPointer		   InputImageConstPointer;
  typedef typename InputImageType::Pointer		   InputImagePointer;
  typedef typename RealImageType::Pointer		   RealImagePointer;
  
  /** Enum of Fast Marching algorithm point types. FarPoints represent far
   * away points; TrialPoints represent points within a narrowband of the
   * propagating front; and KnownPoints represent points which have already
   * been processed. The label NotPoint is given to pixels out of the Final
   * Mask that controls the region where the method is computed.*/
  typedef typename Superclass::LabelType		LabelType;
  
  /** LabelImage typedef support. */
  typedef typename Superclass::LabelPixelType		LabelPixelType;
  typedef typename Superclass::LabelImageType		LabelImageType;
  typedef typename LabelImageType::RegionType		LabelRegionType;
  typedef typename LabelImageType::Pointer		LabelImagePointer;
  typedef typename LabelImageType::OffsetType		OffsetType;

  /** Vector container typedefs */
  typedef typename Superclass::IndexContainerType	IndexContainerType;
  typedef typename Superclass::ThreadIndexContainerType	ThreadIndexContainerType;
  typedef typename Superclass::OffsetContainerType	OffsetContainerType;
  typedef typename Superclass::DoubleContainerType	DoubleContainerType;
  
  /** Required variables to compute values for each stencil */
  typedef typename Superclass::StencilOffsetContainerType	StencilOffsetContainerType;
  typedef typename Superclass::StencilMatrixType		StencilMatrixType;
  typedef typename Superclass::StencilMatrixContainerType	StencilMatrixContainerType;
  typedef typename Superclass::StencilNormsContainerType	StencilNormsContainerType;
  typedef int 		StencilPointerElementType[InputDimension][InputDimension];
  
  /** The image input methods are defined in the base class*/

  
  /** Enum type for choosing the number of stencils, different for 2D and 3D */
  typedef typename Superclass::UseStencilsType		UseStencilsType;

  /** The parameters get/set methods are defined in the base class */

protected:
  MultiStencilEuclideanFastMarching1stOrderImageFilter();
  virtual ~MultiStencilEuclideanFastMarching1stOrderImageFilter() {}
  
  //Overloaded methods
  /** Trial points are stored in a min-heap. This allow efficient access
   * to the trial point with minimum value which is the next grid point
   * the algorithm processes. */
  typedef typename Superclass::NodeType		NodeType;
  typedef typename Superclass::HeapContainer	HeapContainer;
  typedef typename Superclass::NodeComparer	NodeComparer;
  typedef typename Superclass::HeapType		HeapType;
			       
  /** Typenames for helper variables */
  typedef typename Superclass::UpwindVectorType	UpwindVectorType;
  typedef typename Superclass::SignVectorType	SignVectorType;
  
  
  /** Computes the value at the specified index */
  virtual bool ComputeNodeValue( const InputIndexType& idx );
  
  /** CalculaNodosUpwindStencil. Returns true if there is at least one frozen upwind node */
  virtual bool ComputeStencilUpwindNodes( unsigned int stencil,
				  const InputIndexType& idx,
				  UpwindVectorType &upwind,
				  SignVectorType &sign ); 

private:
  MultiStencilEuclideanFastMarching1stOrderImageFilter(const Self&);	// purposely not implemented
  void operator=(const Self&);             // purposely not implemented
  
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMultiStencilEuclideanFastMarching1stOrderImageFilter.txx"
#endif

#endif

