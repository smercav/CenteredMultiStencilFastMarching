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
#ifndef _itkMultiStencilEuclideanFastMarchingImageFilterBase_h
#define _itkMultiStencilEuclideanFastMarchingImageFilterBase_h

#include <functional>
#include <queue>
#include <ctime>
#include <cstring>
#include <string.h>

#include "itkImageToImageFilter.h"
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

namespace itk
{
/** \class MultiStencilEuclideanFastMarchingImageFilterBase
 * \brief Implementation of the configurable multistencil Fast Marching method.
 * To write
 */
template <class TInputImage, class TRealImage, class TOutputImage=TRealImage>
class ITK_EXPORT MultiStencilEuclideanFastMarchingImageFilterBase
		  : public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  /** Convenient typedefs for simplifying declarations. */
  typedef TInputImage                                          InputImageType;
  typedef TRealImage                                           RealImageType;
  typedef TOutputImage                                         OutputImageType;

  /** Standard class typedefs. */
  typedef MultiStencilEuclideanFastMarchingImageFilterBase		Self;
  typedef ImageToImageFilter< InputImageType, OutputImageType> 		Superclass;
  typedef SmartPointer<Self>                                   		Pointer;
  typedef SmartPointer<const Self>                             		ConstPointer;

  /** Template image dimensions */
  static const unsigned int InputDimension	= TInputImage::ImageDimension;
  static const unsigned int RealDimension	= TRealImage::ImageDimension;
  static const unsigned int OutputDimension	= TOutputImage::ImageDimension;
  
  /** Method for creation through the object factory. */
  /*itkNewMacro(Self);*/

  /** Run-time type information (and related methods). */
  itkTypeMacro(MultiStencilEuclideanFastMarchingImageFilter, ImageToImageFilter );

  /** Image typedef support. */
  typedef typename InputImageType::PixelType               InputPixelType;
  typedef typename RealImageType::PixelType                RealPixelType;
  typedef typename OutputImageType::PixelType              OutputPixelType;
  typedef typename NumericTraits<InputPixelType>::RealType InputRealType;
  typedef typename InputImageType::RegionType              InputImageRegionType;
  typedef typename OutputImageType::RegionType             OutputImageRegionType;
  typedef typename OutputImageType::SizeType               OutputSizeType;
  typedef typename InputImageType::SizeType                InputSizeType;
  typedef typename InputImageType::SpacingType		   InputSpacingType;
  typedef typename RealImageType::SizeType                 RealSizeType;
  typedef typename InputImageType::IndexType               InputIndexType;
  typedef typename InputImageType::ConstPointer		   InputImageConstPointer;
  typedef typename InputImageType::Pointer		   InputImagePointer;
  typedef typename RealImageType::Pointer		   RealImagePointer;
  typedef typename OutputImageType::Pointer		   OutputImagePointer;
  
  /** Enum of Fast Marching algorithm point types. FarPoints represent far
   * away points; TrialPoints represent points within a narrowband of the
   * propagating front; and KnownPoints represent points which have already
   * been processed. The label NotPoint is given to pixels out of the Final
   * Mask that controls the region where the method is computed.*/
  typedef enum { FarPoint, KnownPoint, TrialPoint, NotPoint } LabelType;
  
  /** LabelImage typedef support. */
  typedef unsigned int					LabelPixelType;
  typedef Image<LabelPixelType,InputDimension>		LabelImageType;
  typedef typename LabelImageType::RegionType		LabelRegionType;
  typedef typename LabelImageType::Pointer		LabelImagePointer;
  typedef typename LabelImageType::OffsetType		OffsetType;

  /** Vector container typedefs */
//  typedef itk::Array2D<unsigned long>                 ThreadHistogramType;
  typedef std::vector< InputIndexType >			IndexContainerType;
  typedef std::vector< IndexContainerType >		ThreadIndexContainerType;
  typedef std::vector< OffsetType >			OffsetContainerType;
  typedef std::vector< double >				DoubleContainerType;
  
  /** Required variables to compute values for each stencil */
  typedef std::vector< std::vector<OffsetType> >	StencilOffsetContainerType;
  typedef Matrix< double, 
		  InputDimension, 
		  InputDimension >			StencilMatrixType;
  typedef std::vector< StencilMatrixType >		StencilMatrixContainerType;
  typedef std::vector< DoubleContainerType >		StencilNormsContainerType;
  typedef int 		StencilPointerElementType[InputDimension][InputDimension];
  
  /** Image input methods */
  itkSetInputMacro( InitialMaskInput, InputImageType );
  itkGetInputMacro( InitialMaskInput, InputImageType );

  itkSetInputMacro( FinalMaskInput, InputImageType );
  itkGetInputMacro( FinalMaskInput, InputImageType );
  
  itkSetInputMacro( CostInput, RealImageType );
  itkGetInputMacro( CostInput, RealImageType );

  itkSetInputMacro( InitialValuesInput, OutputImageType );
  itkGetInputMacro( InitialValuesInput, OutputImageType );

  /** Filter outputs access methods*/
  itkGetConstReferenceMacro( LabelImage, LabelImagePointer );
  
  virtual const LabelImageType * GetChosenStencilOutput()
  {
    return  const_cast< const LabelImageType * >( m_ChosenStencilImage.GetPointer() );
  }

  
  /** Enum type for choosing the number of stencils. The options are:
   * Use1Stencil,Use4Stencils3D,UseAllStencils */
// #if InputDimension == 2  
//   enum UseStencilsType { OrthoRSquare1, OrthoRSquare2, OrthoRSquare5, OrthoRSquare10, 
// 			 Triangles };
// #else
  enum UseStencilsType { OrthoRSquare1, OrthoRSquare2, OrthoRSquare5, OrthoRSquare6, 
                         OrthoRSquare9, OrthoRSquare10, OrthoRSquare13, 
			 Triangles, Hassouna6Stcl, HassounaExtended };
// #endif
  /** Parameters */
  itkSetMacro( Mask1Value, InputPixelType );
  itkGetConstReferenceMacro( Mask1Value, InputPixelType );
  
  itkSetMacro( Mask2Value, InputPixelType );
  itkGetConstReferenceMacro( Mask2Value, InputPixelType );
  
  itkSetMacro( UseStencils, UseStencilsType );
  itkGetConstReferenceMacro( UseStencils, UseStencilsType );
  
  itkSetMacro( UseFullNhood, bool );
  itkGetConstReferenceMacro( UseFullNhood, bool );

protected:
  MultiStencilEuclideanFastMarchingImageFilterBase();
  virtual ~MultiStencilEuclideanFastMarchingImageFilterBase() {}
  
  //Overloaded methods
  virtual void BeforeThreadedGenerateData();
  virtual void ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, 
#if ITK_VERSION_MAJOR<4
				     int threadId );
#else
				     ThreadIdType threadId );
#endif
  virtual void AfterThreadedGenerateData();
  
  /** Trial points are stored in a min-heap. This allow efficient access
   * to the trial point with minimum value which is the next grid point
   * the algorithm processes. */
  typedef LevelSetNode<OutputPixelType,InputDimension>	NodeType;
  typedef std::vector< NodeType>			HeapContainer;
  typedef std::greater<NodeType>			NodeComparer;
  typedef std::priority_queue< NodeType, 
			       HeapContainer, 
			       NodeComparer >		HeapType;
			       
  /** Typenames for helper variables */
  typedef FixedArray< OutputPixelType, RealDimension >	UpwindVectorType;
  typedef FixedArray< int, RealDimension >		SignVectorType;
  
  /** Typenames for helper classes */
  typedef LinearInterpolateImageFunction<RealImageType>	InterpolatorType;
  typedef typename InterpolatorType::Pointer		InterpolatorPointerType;
  
  typedef typename InterpolatorType::ContinuousIndexType	ContinuousIndexType;
  
  /** Computes the value at the specified index */
  virtual bool ComputeNodeValue( const InputIndexType& idx )=0;
  
  /** Initializes the neighborhood offsets */
  virtual void InitializeNeighborOffsets( bool fullNhood );
  
  /** Computes the stencil coefficients for the scar and the myocardium,
   * returning them by means of pointers */
  virtual void ComputeStencilEquationCoefficients( unsigned int stencil, 
					   const UpwindVectorType &eqUpwind,
					   const SignVectorType& sign,
					   const DoubleContainerType &eqNorms,
					   const double eqCost,
					   double *f0, double *f1, double *f2 );
  
  /** Solve quadratic equation and return the highest root. Returns false if the roots are complex. */
  virtual bool SolveQuadraticEquation( double f0, double f1, double f2, double *sol );
  
  /** Depending on the image dimension and the required number of stencils
   * set the stencil vector offsets, inverse matrices and offset norms*/
  virtual void InitializeStencilElements();
  
  /** CalculaNodosUpwindStencil. Returns true if there is at least one frozen upwind node */
  virtual bool ComputeStencilUpwindNodes( unsigned int stencil,
				  const InputIndexType& idx,
				  UpwindVectorType &upwind,
				  SignVectorType &sign ); 

				
  InputPixelType	m_Mask1Value;
  InputPixelType	m_Mask2Value;
  
  OutputPixelType	m_LargeValue;
  double		m_NumericTol;
  
  InputImageRegionType	m_LargestRegion;
  
  InputSizeType		m_NeighborhoodRadius;
  OffsetContainerType	m_NeighborOffsets;
  DoubleContainerType	m_NeighborNorms;
  bool			m_UseFullNhood;
  
  LabelImagePointer	m_LabelImage;
  InputImagePointer	m_FinalMaskImage;
  LabelImagePointer	m_ChosenStencilImage;
  RealImagePointer	m_CostImage;	//Auxiliary pointer to improve performance
  OutputImagePointer	m_InitialValuesImage;
  
  
  InterpolatorPointerType	m_Interpolator;
  
  ThreadIndexContainerType	m_NarrowBandIndexContainer;
  HeapType			m_TrialHeap;
  
  StencilOffsetContainerType	m_StencilOffsets;
  StencilMatrixContainerType	m_StencilInvMatrices;
  StencilNormsContainerType	m_StencilOffsetNorms;
  
  UseStencilsType 		m_UseStencils;
  unsigned int			m_NumberOfStencils;
  StencilPointerElementType const*	m_StencilCoefficients;

  //These are the general stencil sets that may be used
  static const int		m_StencilCoefficients2DOrthogonal[6][2][2];
  static const int		m_StencilCoefficients3DOrthogonal[38][3][3];
  static const int		m_StencilCoefficients3DHassouna[10][3][3];
  static const int		m_StencilCoefficients2DTriangles[4][2][2];
  static const int		m_StencilCoefficients3DTriangles[24][3][3];
  
  clock_t			m_RunningStartTime;
  double			m_RunningSeconds;
  
private:
  // Constructor and operator= purposely not implemented
  MultiStencilEuclideanFastMarchingImageFilterBase(const Self&);
  void operator=(const Self&);

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMultiStencilEuclideanFastMarchingImageFilterBase.txx"
#endif

#endif

