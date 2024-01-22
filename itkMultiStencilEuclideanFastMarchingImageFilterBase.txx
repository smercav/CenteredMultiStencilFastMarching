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
#ifndef _itk_MultiStencilEuclideanFastMarchingImageFilterBase_txx
#define _itk_MultiStencilEuclideanFastMarchingImageFilterBase_txx
#include "itkMultiStencilEuclideanFastMarchingImageFilterBase.h"

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


template <class TInputImage, class TRealImage, class TOutputImage>
MultiStencilEuclideanFastMarchingImageFilterBase<TInputImage, TRealImage, TOutputImage>
::MultiStencilEuclideanFastMarchingImageFilterBase() 
{
  //Change the name of the primary input so that the [SG]etScarInput() methods work
  typename Superclass::DataObjectIdentifierType currPrimaryInputName( "InitialMaskInput" );
  this->Superclass::SetPrimaryInputName( currPrimaryInputName );
  this->ProcessObject::AddRequiredInputName( std::string("CostInput") );
    
//   this->SetNumberOfRequiredOutputs( 2 );
//   this->SetNthOutput( 0, this->MakeOutput(0) );  
//   this->SetNthOutput( 1, this->MakeOutput(1) );
  
  //Set parameters
  m_Mask1Value = NumericTraits<InputPixelType>::max();
  m_Mask2Value = m_Mask1Value;
  m_LargeValue = NumericTraits<OutputPixelType>::max();
  m_NeighborhoodRadius.Fill( 1 );
  
  m_Interpolator = InterpolatorType::New();
  
  m_UseStencils = OrthoRSquare2;
  m_UseFullNhood = false;
  
  m_NumericTol = 1.0e-10f;
  
  //this->SetNumberOfThreads( 1 );
}

template< class TInputImage, class TRealImage, class TOutputImage>
void 
MultiStencilEuclideanFastMarchingImageFilterBase< TInputImage, TRealImage, TOutputImage>
::BeforeThreadedGenerateData(  ) 
{
    // Do this first:
    this->Superclass::BeforeThreadedGenerateData();
    
    //Get the start time
    m_RunningStartTime = clock();

    //Check that the inputs are ok
    if ( (InputDimension!=RealDimension)||(RealDimension!=OutputDimension) )
    {
      itkExceptionMacro(<< "Different input dimensions.");
    }
    //else if ( (InputDimension<2) || (InputDimension>3) )
    else if ((InputDimension<2)||(InputDimension>3))
    {
      itkExceptionMacro(<< "This filter is only supported for 2D and 3D images." );
//       itkExceptionMacro(<< "This filter is only supported for 3D images." );
    }
    else
    {
      InputSizeType sz0 = this->GetInitialMaskInput()->GetLargestPossibleRegion().GetSize();
      RealSizeType sz1 = this->GetCostInput()->GetLargestPossibleRegion().GetSize();
      InputSizeType sz2 = (this->GetFinalMaskInput())? this->GetFinalMaskInput()->GetLargestPossibleRegion().GetSize() : sz0;
      OutputSizeType sz3 = (this->GetInitialValuesInput())? this->GetInitialValuesInput()->GetLargestPossibleRegion().GetSize() : sz0;
      
      for (unsigned int i=0; i<InputDimension; i++)
      {
	if ( (sz0[i]!=sz1[i]) || (sz1[i]!=sz2[i]) || (sz1[i]!=sz3[i]) )
	{
	  itkExceptionMacro(<< "There are inputs with different sizes.");
	}
      }
    }
    
    //If the final mask image is not set, an all-ones image is created.
    if (!(this->GetFinalMaskInput()))
    {
      m_FinalMaskImage = InputImageType::New();
      m_FinalMaskImage->SetDirection( this->GetInitialMaskInput()->GetDirection() );
      m_FinalMaskImage->SetOrigin( this->GetInitialMaskInput()->GetOrigin() );
      m_FinalMaskImage->SetSpacing( this->GetInitialMaskInput()->GetSpacing() );
      m_FinalMaskImage->SetRegions( this->GetInitialMaskInput()->GetLargestPossibleRegion() );
      m_FinalMaskImage->Allocate();
      m_FinalMaskImage->FillBuffer( m_Mask2Value );
    }
    else
    {
      m_FinalMaskImage = const_cast< InputImageType * >( this->GetFinalMaskInput() );
    }
    
    //If the init values image is not set, an all-zeros image is created.
    if (!(this->GetInitialValuesInput()))
    {
      m_InitialValuesImage = OutputImageType::New();
      m_InitialValuesImage->SetDirection( this->GetCostInput()->GetDirection() );
      m_InitialValuesImage->SetOrigin( this->GetCostInput()->GetOrigin() );
      m_InitialValuesImage->SetSpacing( this->GetCostInput()->GetSpacing() );
      m_InitialValuesImage->SetRegions( this->GetCostInput()->GetLargestPossibleRegion() );
      m_InitialValuesImage->Allocate();
      m_InitialValuesImage->FillBuffer( NumericTraits<OutputPixelType>::Zero );
    }
    else
    {
      m_InitialValuesImage = const_cast< OutputImageType * >( this->GetInitialValuesInput() );
    }
    
    //Build the thread container for the narrow band nodes and reset it
    unsigned int k = this->GetNumberOfThreads();
    m_NarrowBandIndexContainer.resize( k );
    for (unsigned int i=0; i<k; i++)
    {
      m_NarrowBandIndexContainer.at(i).clear();
    }
    
    //Create the neighbor offset table. 
    this->InitializeNeighborOffsets( m_UseFullNhood );
    
    //compute the neighbor norms, spacing included
    InputSpacingType spacing = this->GetInitialMaskInput()->GetSpacing();
    m_NeighborNorms.clear();
    typename OffsetContainerType::const_iterator it;
    
    for (it=m_NeighborOffsets.begin(); it!=m_NeighborOffsets.end(); it++)
    {
      double norm = NumericTraits<double>::Zero;
      for (unsigned int i=0; i<InputDimension; i++)
      {
	norm += std::pow( (*it)[i]*spacing[i], 2 );
      }
      m_NeighborNorms.push_back( std::sqrt( norm ) );
    }
    
    //Compute the stencil vectors, offsets, etc
    this->InitializeStencilElements();
    
    //Allocate the label image
    InputImageConstPointer input  = this->GetInput();
    m_LabelImage = LabelImageType::New();
    m_LabelImage->SetOrigin(  input->GetOrigin() );
    m_LabelImage->SetSpacing( input->GetSpacing() );
    m_LabelImage->SetRegions( input->GetLargestPossibleRegion() );
    m_LabelImage->SetDirection( input->GetDirection() );
    m_LabelImage->Allocate();
    
    //Allocate the chosen stencil image
    m_ChosenStencilImage = LabelImageType::New();
    m_ChosenStencilImage->SetOrigin(  input->GetOrigin() );
    m_ChosenStencilImage->SetSpacing( input->GetSpacing() );
    m_ChosenStencilImage->SetRegions( input->GetLargestPossibleRegion() );
    m_ChosenStencilImage->SetDirection( input->GetDirection() );
    m_ChosenStencilImage->Allocate();
    m_ChosenStencilImage->FillBuffer( m_NumberOfStencils );
    
    //Get the largest image region
    m_LargestRegion = input->GetLargestPossibleRegion();
    
    // Initialize the InterpolatorPointerType
    m_Interpolator->SetInputImage( this->GetCostInput() );
    
    // To avoid calling GetCostInput so many times, we store the pointer in a class member
    m_CostImage = const_cast<RealImageType *>( this->GetCostInput() );

    itkDebugMacro("BeforeThreadedGenerateData complete");
}

template< class TInputImage, class TRealImage, class TOutputImage>
void 
MultiStencilEuclideanFastMarchingImageFilterBase< TInputImage, TRealImage, TOutputImage>
#if ITK_VERSION_MAJOR<4
::ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, int threadId ) 
#else
::ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, ThreadIdType threadId ) 
#endif
{
  //Here, the image is iterated to find the pixels with 
  //at least one known neighbor and add it to the narrow band
  //The iteration is also exploited to initialize 
  //the label image and the outputs.
  itkDebugMacro( "Beginning ThreadedGenerateData " << threadId );
  
  //Iterator typedefs 
  typedef ImageRegionIterator<OutputImageType>		OutputRegionIterator;
  typedef ImageRegionIterator<LabelImageType>		LabelRegionIterator;
  typedef ConstNeighborhoodIterator<InputImageType>	InputNeighborhoodIterator;
  typedef typename InputNeighborhoodIterator
				::NeighborhoodType	NeighborhoodType;
  typedef NeighborhoodAlgorithm
      ::ImageBoundaryFacesCalculator<OutputImageType>	FaceCalculatorType;
  typedef typename FaceCalculatorType::FaceListType	FaceListType;
  

  //Scan the image, initialize the label and the output images, 
  //and store the narrow band indices, iterating on each of 
  //the divided regions
  FaceCalculatorType faceCalculator;
  FaceListType faceList = faceCalculator(this->GetOutput(0), 
					 outputRegionForThread, m_NeighborhoodRadius );
  
  typename FaceListType::iterator fit;
  for (fit=faceList.begin(); fit!=faceList.end(); ++fit)
  {
    //Define the iterators using this regions
    InputNeighborhoodIterator imIt( m_NeighborhoodRadius, this->GetInitialMaskInput(), *fit );
    InputNeighborhoodIterator fmIt( m_NeighborhoodRadius, this->m_FinalMaskImage, *fit );
    LabelRegionIterator lIt( m_LabelImage, *fit );
    OutputRegionIterator tmIt( this->GetOutput(), *fit );
    OutputRegionIterator ivIt( m_InitialValuesImage, *fit );
    
    imIt.GoToBegin(); fmIt.GoToBegin(); lIt.GoToBegin(); tmIt.GoToBegin(); ivIt.GoToBegin();
    while (!imIt.IsAtEnd())
    {
      //Check the initial mask central pixel, depending on its value
      //initialize the labels, outputs, and narrow band
      if ( imIt.GetCenterPixel() == m_Mask1Value )
      {
	//Set the values
	lIt.Set( KnownPoint );
	tmIt.Set( ivIt.Value() );
      }
      else
      {
	//Set the values
	if (fmIt.GetCenterPixel()==m_Mask2Value)
	{
	  tmIt.Set( m_LargeValue );
	  lIt.Set( FarPoint );
	
	  //Check if it is a narrow band pixel. If so,
	  //add its index to the container
	  //InputIndexType ctrIdx = imIt.GetIndex(); //For debugging
	  
	  for (unsigned int i=0; i<(unsigned int)imIt.Size(); i++)
	  {
	    if ( (imIt.GetPixel(i)==m_Mask1Value) && (fmIt.GetPixel(i)==m_Mask2Value) )
	    {
	      m_NarrowBandIndexContainer.at(threadId).push_back( imIt.GetIndex() );
	      break;
	    }
	  } //End for neighborhood
	  
	}
	else
	{
	  //The point is outside the mask where we wish to compute
	  //the algorithm. It is initialized as Zero and Not a point
	  tmIt.Set( NumericTraits<OutputPixelType>::Zero );
	  lIt.Set( NotPoint );
	}
	
      }
      
      //Increase iterators
      ++imIt; ++fmIt; ++lIt; ++tmIt; ++ivIt;
    }
  } //End faceList FOR

  itkDebugMacro( "ThreadedGenerateData " << threadId << " complete." );
  
}

template< class TInputImage, class TRealImage, class TOutputImage>
void 
MultiStencilEuclideanFastMarchingImageFilterBase< TInputImage, TRealImage, TOutputImage>
::AfterThreadedGenerateData() 
{
  
  itkDebugMacro( "Beginning AfterThreadedGenerateData" );

  typedef ConstNeighborhoodIterator<LabelImageType>	  NeighborhoodIterator;
  typedef typename NeighborhoodIterator::NeighborhoodType NeighborhoodType;

  //Convenient image pointers
  typename OutputImageType::Pointer Tm = this->GetOutput();
  
  //Create a NeighborhoodIterator that will be used to get 
  //pixel neighborhoods at random locations. It may not be
  //computationally efficient.
  //LabelRegionType m_LargestRegion = m_LabelImage->GetLargestPossibleRegion();
  NeighborhoodIterator labIt( m_NeighborhoodRadius, m_LabelImage, m_LargestRegion );
  labIt.NeedToUseBoundaryConditionOff();
  
  //Build the priority queue from the narrow band index container
  //populated in ThreadedGenerateData. Compute its values too.
  // Make sure the heap is empty (from itk::FastMarchingImageFilter)
  m_TrialHeap = HeapType();
  while ( !m_TrialHeap.empty() )
  {
    m_TrialHeap.pop();
  }
  
  
  for ( unsigned int k=0; k<m_NarrowBandIndexContainer.size(); k++ )
  {
    //Traverse the index vector
    typename IndexContainerType::const_iterator it;
    for (it  = m_NarrowBandIndexContainer.at(k).begin();
	 it != m_NarrowBandIndexContainer.at(k).end(); it++)
    {
      
      //Compute its value. There may be instances where no valid solution is found.
      //If that is the case, that node is not incorporated to the heap.
      if (this->ComputeNodeValue( *it ))
      {

	//Change the label to active
	m_LabelImage->SetPixel( *it, TrialPoint );
	
	//Insert the node into the heap 
	OutputPixelType fvalue = Tm->GetPixel( *it );
	NodeType newNode; 
	newNode.SetValue(fvalue); newNode.SetIndex(*it);
	m_TrialHeap.push( newNode );
      }
    }
  }
  itkDebugMacro( "Narrow band set up." );
  
  itkDebugMacro( "Starting region growing step..." );
  //While there are elements in the priority queue, perform the
  //main loop of the Streamline Fast Marching algorithm. However, now 
  //it's not a good idea to use a ConstNeighborhoodIterator, because
  //the offsets and indices of the neighbors are needed.
  unsigned int neighSize = m_NeighborOffsets.size();
  OutputPixelType       currNodeValue;
  
  while (!m_TrialHeap.empty())
  {
    //Get the node to be processed next,
    //and take it out of the heap.
    //It should be taken into account that there may be
    //repeated nodes in the heap.
    
    InputIndexType currIdx  = m_TrialHeap.top().GetIndex();
    OutputPixelType currValue = m_TrialHeap.top().GetValue();
    m_TrialHeap.pop();
    
    // does this node contain the current value ?
    currNodeValue = Tm->GetPixel( currIdx );
    
    itkDebugMacro("Considered IDX, value and nodevalue: " << currIdx << " ; " << currValue << ", " << currNodeValue );

    if ( (Math::ExactlyEquals(currValue, currNodeValue)) && (m_LabelImage->GetPixel(currIdx)==TrialPoint) )
    {
      itkDebugMacro("Processing IDX: " << currIdx );
      
      if ((currIdx[0]==157)&&(currIdx[1]==137)&&(currIdx[2]==4))
	std::cout << "Nodo a depurar alcanzado" << std::endl;
      
      //Label it known 
      m_LabelImage->SetPixel( currIdx, KnownPoint );
      
      //Get its neighbors and process them if they are labeled
      //Active or Far within the mask
      for (unsigned int i=0; i<neighSize; i++)
      {
	InputIndexType neighIdx = currIdx + m_NeighborOffsets[i];
	
	if (m_LargestRegion.IsInside(neighIdx)) 
	{
	  LabelPixelType label = m_LabelImage->GetPixel(neighIdx);
	  if (label==FarPoint)
	  {
	    //Compute its value
	    if (this->ComputeNodeValue( neighIdx ))
	    {
	      //Change its label to TrialPoint
	      m_LabelImage->SetPixel( neighIdx, TrialPoint );
	      
	      //Insert the node into the heap 
	      OutputPixelType fvalue = Tm->GetPixel( neighIdx );
	      NodeType newNode; 
	      newNode.SetValue(fvalue); newNode.SetIndex(neighIdx);
	      m_TrialHeap.push( newNode );
	    }
	  }
	  else if (label==TrialPoint)
	  {
	    OutputPixelType oldValue = Tm->GetPixel( neighIdx );
	    
	    //Compute its value
	    this->ComputeNodeValue( neighIdx );
	    
	    OutputPixelType fvalue = Tm->GetPixel( neighIdx );
	    
	    if ( oldValue > fvalue )
	    {
	      //Insert the node into the heap 
	      NodeType newNode; 
	      newNode.SetValue(fvalue); newNode.SetIndex(neighIdx);
	      m_TrialHeap.push( newNode );
	    }
	    else if (oldValue < fvalue)
	    {
	      //The previous value was lower. Undo the change...
	      Tm->SetPixel( neighIdx, oldValue );
	    }
	  }
	}
      }
    
    }
  }

  clock_t endingTime = clock();
  m_RunningSeconds = static_cast<double>( (endingTime-m_RunningStartTime)/CLOCKS_PER_SEC );
  std::cout << "Running seconds = " << m_RunningSeconds << std::endl;
  
  //Do this last:
  this->Superclass::AfterThreadedGenerateData();
}





template< class TInputImage, class TRealImage, class TOutputImage>
void 
MultiStencilEuclideanFastMarchingImageFilterBase< TInputImage, TRealImage, TOutputImage>
::InitializeNeighborOffsets( bool fullNhood )
{

    //Create the neighbor offset table. A const zero offset is needed
    //to be able to create negative offsets by substraction.
    OffsetType offZero; offZero.Fill( 0 ); 
    m_NeighborOffsets.clear();
    for (unsigned int i=0; i<InputDimension; i++)
    {
      OffsetContainerType prevOffsets( m_NeighborOffsets );

      OffsetType offDim = OffsetType::GetBasisOffset( i );
      m_NeighborOffsets.push_back( offDim );
      m_NeighborOffsets.push_back( offZero-offDim );
      
      if (fullNhood)
      {
	for (unsigned int j=0; j<prevOffsets.size(); j++)
	{
	  m_NeighborOffsets.push_back( prevOffsets[j]+offDim );
	  m_NeighborOffsets.push_back( prevOffsets[j]-offDim );
	}
      }
    }
  
}


template< class TInputImage, class TRealImage, class TOutputImage>
void 
MultiStencilEuclideanFastMarchingImageFilterBase< TInputImage, TRealImage, TOutputImage>
::InitializeStencilElements()
{
  typedef typename StencilMatrixType::InternalMatrixType StencilInternalMatrixType;
  
  unsigned int numVecXSten = InputDimension;
  
  //Determine how many stencils we need to use, 
  //unsigned int numStencils; //in class variable m_NumberOfStencils
if (InputDimension==2)
{
  if ((m_UseStencils==OrthoRSquare1)||(m_UseStencils==OrthoRSquare2)
    ||(m_UseStencils==OrthoRSquare5)||(m_UseStencils==OrthoRSquare10))
  {
    switch (m_UseStencils)
    {
      case OrthoRSquare1:
	m_NumberOfStencils = 1;
	break;
      case OrthoRSquare2:
	m_NumberOfStencils = 2;
	break;
      case OrthoRSquare5:
	m_NumberOfStencils = 4;
	break;
      case OrthoRSquare10:
	m_NumberOfStencils = 6;
	break;
    }
    
    m_StencilOffsets.resize(m_NumberOfStencils);
    for (unsigned int stencil=0; stencil<m_NumberOfStencils; stencil++)
    {
      //Initialize variables
      m_StencilOffsets[stencil].clear();
      for (unsigned int i=0; i<numVecXSten; i++)
      {
	OffsetType offset;
	
	for (unsigned int j=0; j<InputDimension; j++)
	{
	  offset[j] = m_StencilCoefficients2DOrthogonal[stencil][i][j];
	}
	m_StencilOffsets[stencil].push_back( offset );
      }
    }
  }
  else if (m_UseStencils==Triangles)
  {
    m_NumberOfStencils = 4;
  
    m_StencilOffsets.resize(m_NumberOfStencils);
    for (unsigned int stencil=0; stencil<m_NumberOfStencils; stencil++)
    {
      //Initialize variables
      m_StencilOffsets[stencil].clear();
      for (unsigned int i=0; i<numVecXSten; i++)
      {
	OffsetType offset;
	
	for (unsigned int j=0; j<InputDimension; j++)
	{
	  offset[j] = m_StencilCoefficients2DTriangles[stencil][i][j];
	}
	m_StencilOffsets[stencil].push_back( offset );
      }
    }
  }
  else
  {
    itkExceptionMacro( << "Unknown stencil set." );
  }
}
else if (InputDimension==3)
{
  if ((m_UseStencils==OrthoRSquare1)||(m_UseStencils==OrthoRSquare2)||(m_UseStencils==OrthoRSquare5)
    ||(m_UseStencils==OrthoRSquare6)||(m_UseStencils==OrthoRSquare9)
    ||(m_UseStencils==OrthoRSquare10)||(m_UseStencils==OrthoRSquare13))
  {
  
    switch (m_UseStencils)
    {
      case OrthoRSquare1:	
	m_NumberOfStencils = 1;
	break;
      case OrthoRSquare2:
	m_NumberOfStencils = 4;
	break;
      case OrthoRSquare5:
	m_NumberOfStencils = 10;
	break;
      case OrthoRSquare6:
	m_NumberOfStencils = 22;
	break;
      case OrthoRSquare9:
	m_NumberOfStencils = 26;
	break;
      case OrthoRSquare10:
	m_NumberOfStencils = 32;
	break;
      case OrthoRSquare13:
	m_NumberOfStencils = 38;
	break;
    }
    
    m_StencilOffsets.resize(m_NumberOfStencils);
    for (unsigned int stencil=0; stencil<m_NumberOfStencils; stencil++)
    {
      //Initialize variables
      m_StencilOffsets[stencil].clear();
      for (unsigned int i=0; i<numVecXSten; i++)
      {
	OffsetType offset;
	
	for (unsigned int j=0; j<InputDimension; j++)
	{
	  offset[j] = m_StencilCoefficients3DOrthogonal[stencil][i][j];
	}
	m_StencilOffsets[stencil].push_back( offset );
      }
    }
  }
  else if ((m_UseStencils==Hassouna6Stcl)||(m_UseStencils==HassounaExtended))
  {
  
    switch (m_UseStencils)
    {
      case Hassouna6Stcl:
	m_NumberOfStencils = 6;
	break;
      case HassounaExtended:
	m_NumberOfStencils = 10;
	break;
    }
    
    m_StencilOffsets.resize(m_NumberOfStencils);
    for (unsigned int stencil=0; stencil<m_NumberOfStencils; stencil++)
    {
      //Initialize variables
      m_StencilOffsets[stencil].clear();
      for (unsigned int i=0; i<numVecXSten; i++)
      {
	OffsetType offset;
	
	for (unsigned int j=0; j<InputDimension; j++)
	{
	  offset[j] = m_StencilCoefficients3DHassouna[stencil][i][j];
	}
	m_StencilOffsets[stencil].push_back( offset );
      }
    }
	
  }
  else if (m_UseStencils==Triangles)
  {
    m_NumberOfStencils = 24;
  
    m_StencilOffsets.resize(m_NumberOfStencils);
    for (unsigned int stencil=0; stencil<m_NumberOfStencils; stencil++)
    {
      //Initialize variables
      m_StencilOffsets[stencil].clear();
      for (unsigned int i=0; i<numVecXSten; i++)
      {
	OffsetType offset;
	
	for (unsigned int j=0; j<InputDimension; j++)
	{
	  offset[j] = m_StencilCoefficients3DTriangles[stencil][i][j];
	}
	m_StencilOffsets[stencil].push_back( offset );
      }
    }
  }
  else
  {
    itkExceptionMacro( << "Unknown stencil set." );
  }
}
else
{
  itkExceptionMacro(<< "This filter is only supported for 2D and 3D images." );
}
  
  //Initialize vectors
  m_StencilOffsets.resize(m_NumberOfStencils);
  m_StencilInvMatrices.resize(m_NumberOfStencils);
  m_StencilOffsetNorms.resize(m_NumberOfStencils);
  
  //Get the image spacing
  InputSpacingType spacing = this->GetInitialMaskInput()->GetSpacing();
  
  for (unsigned int stencil=0; stencil<m_NumberOfStencils; stencil++)
  {
    //Initialize variables
    m_StencilOffsets[stencil].clear();
    m_StencilOffsetNorms[stencil].clear();
    StencilMatrixType auxMatrix;
    
    for (unsigned int i=0; i<numVecXSten; i++)
    {
      //Offsets
      OffsetType offset = m_StencilOffsets[stencil][i];
      double vecNorm = NumericTraits<double>::Zero;
      for (unsigned int j=0; j<InputDimension; j++)
      {
	vecNorm += std::pow( spacing[j]*offset[j], 2 );
      }
      vecNorm = std::sqrt( vecNorm );
      
      m_StencilOffsetNorms[stencil].push_back( vecNorm );
      
      //Rellenamos la fila de la matriz
      for (unsigned int j=0; j<InputDimension; j++)
      {
	auxMatrix[i][j] = (spacing[j]*offset[j]) / vecNorm;
      }
    }
    
    // Multiply by its transpose
    auxMatrix *= auxMatrix.GetTranspose();
    
    //Store the matrix inverse
    StencilInternalMatrixType auxInverse = auxMatrix.GetInverse();
    
    //Copy the inverse into a StencilMatrixType object
    StencilMatrixType inverse;
    for (unsigned int r=0; r<InputDimension; r++)
    {
      for (unsigned int c=0; c<InputDimension; c++)
      { 
	inverse(r,c) = auxInverse.get(r,c);
      }
    }
    
    m_StencilInvMatrices[stencil] = auxInverse; //inverse;
//     std::cout << "Inverse matrix stencil " << stencil << " =" << std::endl;
//     std::cout << m_StencilInvMatrices[stencil] << std::endl;
  }
}


template< class TInputImage, class TRealImage, class TOutputImage>
void 
MultiStencilEuclideanFastMarchingImageFilterBase< TInputImage, TRealImage, TOutputImage>
::ComputeStencilEquationCoefficients( unsigned int stencil, 
				      const UpwindVectorType &eqUpwind,
				      const SignVectorType& sign,
				      const DoubleContainerType &eqNorms,
				      const double eqCost,
				      double *f0, double *f1, double *f2 )
{
  // WARNING: IMPLEMENTADA, SIN DEPURAR
  //Compute the stencil equation coefficients
  *f0 = NumericTraits<double>::Zero;
  *f1 = NumericTraits<double>::Zero;
  *f2 = NumericTraits<double>::Zero;
  
  double Wdiag, Wcross;
  
  for (unsigned int i=0; i<InputDimension; i++)
  {
    if (sign[i]==0)
      continue;
    
    Wdiag = m_StencilInvMatrices[stencil][i][i] / (eqNorms[i]*eqNorms[i]);
    
    *f2 += Wdiag;
    *f1 -= Wdiag * 2*eqUpwind[i];
    *f0 += Wdiag * eqUpwind[i]*eqUpwind[i];
    
    for (unsigned int j=i+1; j<InputDimension; j++)
    {
      if ((sign[j]==0)||(m_StencilInvMatrices[stencil][i][j]==0))
	continue;
      
      Wcross = 2*m_StencilInvMatrices[stencil][i][j] *(sign[i]*sign[j]) / (eqNorms[i]*eqNorms[j]);
      
      *f2 += Wcross;
      *f1 -= Wcross * (eqUpwind[i]+eqUpwind[j]);
      *f0 += Wcross * (eqUpwind[i]*eqUpwind[j]);

    }
  }
  
  //Remove the cost^2 from f0
  *f0 -= eqCost*eqCost;

} //End ComputeStencilEquationCoefficients




template< class TInputImage, class TRealImage, class TOutputImage>
bool 
MultiStencilEuclideanFastMarchingImageFilterBase< TInputImage, TRealImage, TOutputImage>
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
      candval[i] = m_LargeValue;
    
    //Check if the upwind node is within the image and frozen
    fwIdx = idx + m_StencilOffsets[stcl][dim];
    if ((m_LargestRegion.IsInside(fwIdx)) && (m_LabelImage->GetPixel(fwIdx)==KnownPoint))
    {
      candval[0] = output->GetPixel(fwIdx);
    }
    
    //Check if the downwind node is within the image and frozen
    bwIdx = idx - m_StencilOffsets[stcl][dim];
    if ((m_LargestRegion.IsInside(bwIdx)) && (m_LabelImage->GetPixel(bwIdx)==KnownPoint))
    {
      candval[1] = output->GetPixel(bwIdx);
    }
    
    //Depending on candval, assign upwind[dim] and sign[dim]
    if ( (candval[0]==m_LargeValue) && (candval[1]==m_LargeValue) )
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
bool 
MultiStencilEuclideanFastMarchingImageFilterBase< TInputImage, TRealImage, TOutputImage>
::SolveQuadraticEquation( double f0, double f1, double f2, double *sol )
{
  double discri = f1*f1 - 4*f2*f0;
  
  if (discri>0)
  {
    double sol1, sol2;
    
    double sqrDiscri = std::sqrt( discri );
    sol1 = (-f1 + sqrDiscri) / (2*f2);
    sol2 = (-f1 - sqrDiscri) / (2*f2);
    
    *sol = (sol1>sol2)? sol1 : sol2;
    return true;
  }
  else if (discri>-m_NumericTol) //Tolerate almost zero discriminants
  {
    *sol = - f1 / (2*f2);
    return true;
  }
  else
  {
    return false;
  }
    
}



template <class TInputImage, class TRealImage, class TOutputImage>
const int
MultiStencilEuclideanFastMarchingImageFilterBase<TInputImage, TRealImage, TOutputImage>
::m_StencilCoefficients2DOrthogonal[][2][2] = { {{1,0},{0,1}}, {{1,1},{-1,1}},
  {{2,1},{1,-2}},{{2,-1},{1,2}},{{3,1},{1,-3}},{{3,-1},{1,3}}};

template <class TInputImage, class TRealImage, class TOutputImage>
const int
MultiStencilEuclideanFastMarchingImageFilterBase<TInputImage, TRealImage, TOutputImage>
::m_StencilCoefficients3DOrthogonal[][3][3] = {
  { {1,0,0}, {0,1,0}, {0,0,1} }, { {0,0,1}, {1,1,0}, {1,-1,0} },
  { {1,0,0}, {0,1,1}, {0,-1,1}}, { {0,1,0}, {1,0,1}, {-1,0,1} },
  { {0,0,1}, {2,1,0}, {1,-2,0} },{ {0,0,1}, {2,-1,0}, {1,2,0} },
  { {1,0,0}, {0,1,2}, {0,-2,1} },{ {1,0,0}, {0,-1,2},{0,2,1} },
  {{0,1,0},{1,0,2},{-2,0,1}},{{0,1,0},{-1,0,2},{2,0,1}},
  {{1,0,1},{-1,1,1},{-1,-2,1}},{{1,0,1},{-1,-1,1},{-1,2,1}},
  {{0,1,1},{1,-1,1},{-2,-1,1}},{{0,1,1},{-1,-1,1},{2,-1,1}},
  {{0,-1,1},{1,1,1},{-2,1,1}},{{0,-1,1},{-1,1,1},{2,1,1}},
  {{-1,0,1},{1,1,1},{1,-2,1}},{{-1,0,1},{1,-1,1},{1,2,1}},
  {{1,1,0},{1,-1,1},{-1,1,2}},{{1,1,0},{-1,1,1},{1,-1,2}},
  {{1,-1,0},{1,1,1},{-1,-1,2}},{{1,-1,0},{-1,-1,1},{1,1,2}},
  {{2,1,2},{-1,-2,2},{-2,2,1}},{{2,-1,2},{-1,2,2},{-2,-2,1}},
  {{1,2,2},{-2,-1,2},{2,-2,1}},{{1,-2,2},{-2,1,2},{2,2,1}},
  {{0,0,1},{3,1,0},{1,-3,0}},{{0,0,1},{3,-1,0},{1,3,0}},
  {{1,0,0},{0,1,3},{0,-3,1}},{{1,0,0},{0,-1,3},{0,3,1}},
  {{0,1,0},{1,0,3},{-3,0,1}},{{0,1,0},{-1,0,3},{3,0,1}},
  {{0,0,1},{3,2,0},{2,-3,0}},{{0,0,1},{3,-2,0},{2,3,0}},
  {{1,0,0},{0,2,3},{0,-3,2}},{{1,0,0},{0,-2,3},{0,3,2}},
  {{0,1,0},{2,0,3},{-3,0,2}},{{0,1,0},{-2,0,3},{3,0,2}} };

  
template <class TInputImage, class TRealImage, class TOutputImage>
const int
MultiStencilEuclideanFastMarchingImageFilterBase<TInputImage, TRealImage, TOutputImage>
::m_StencilCoefficients3DHassouna[][3][3] = {
  { {1,0,0},{0,1,0},{0,0,1}  }, { {1,0,0},{0,1,-1},{0,1,1} },
  { {1,0,-1},{0,1,0},{1,0,1} }, { {1,1,0},{-1,1,0},{0,0,1} },
  {{1,0,1},{-1,1,1},{1,1,-1}},  {{1,0,-1},{1,1,1},{-1,1,-1}},
  {{1,1,0},{1,-1,1},{-1,1,1}},  {{-1,1,0},{1,1,1},{-1,-1,1}},
  {{0,1,1},{1,1,-1},{1,-1,1}},  {{0,-1,1},{1,1,1},{1,-1,-1}} };

  
template <class TInputImage, class TRealImage, class TOutputImage>
const int
MultiStencilEuclideanFastMarchingImageFilterBase<TInputImage, TRealImage, TOutputImage>
::m_StencilCoefficients2DTriangles[][2][2] = { {{1,0},{1,1}},{{1,1},{0,1}},{{0,1},{-1,1}},{{-1,1},{-1,0}} }; 


template <class TInputImage, class TRealImage, class TOutputImage>
const int
MultiStencilEuclideanFastMarchingImageFilterBase<TInputImage, TRealImage, TOutputImage>
::m_StencilCoefficients3DTriangles[][3][3] = {
{{1,0,0},{1,1,0},{1,1,1}},{{1,0,0},{1,1,1},{1,0,1}},{{1,0,0},{1,0,1},{1,-1,1}},{{1,0,0},{1,-1,1},{1,-1,0}},
{{1,0,0},{1,-1,0},{1,-1,-1}},{{1,0,0},{1,-1,-1},{1,0,-1}},{{1,0,0},{1,0,-1},{1,1,-1}},{{1,0,0},{1,1,-1},{1,1,0}},
{{0,1,0},{0,1,1},{1,1,1}},{{0,1,0},{1,1,1},{1,1,0}},{{0,1,0},{1,1,0},{1,1,-1}},{{0,1,0},{1,1,-1},{0,1,-1}},
{{0,1,0},{0,1,-1},{-1,1,-1}},{{0,1,0},{-1,1,-1},{-1,1,0}},{{0,1,0},{-1,1,0},{-1,1,1}},{{0,1,0},{-1,1,1},{0,1,1}},
{{0,0,1},{1,0,1},{1,1,1}},{{0,0,1},{1,1,1},{0,1,1}},{{0,0,1},{0,1,1},{-1,1,1}},{{0,0,1},{-1,1,1},{-1,0,1}},
{{0,0,1},{-1,0,1},{-1,-1,1}},{{0,0,1},{-1,-1,1},{0,-1,1}},{{0,0,1},{0,-1,1},{1,-1,1}},{{0,0,1},{1,-1,1},{1,0,1}} };

} // end namespace itk

#endif //#ifndef





