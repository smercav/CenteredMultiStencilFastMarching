/*=========================================================================

 =========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif


#include "itkPluginUtilities.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkNumericTraits.h"

#include "itkMultiStencilEuclideanFastMarchingImageFilterBase.h"
#include "itkMultiStencilEuclideanFastMarching1stOrderImageFilter.h"
#include "itkMultiStencilEuclideanFastMarching2ndOrderImageFilter.h"
#include "itkMultiStencilEuclideanFastMarchingCenteredImageFilter.h"

#include <boost/program_options.hpp>
#include <iostream>
#include <string>

namespace po = boost::program_options;


template<class TMask, class TFeature>  
int DoIt( po::variables_map &vm, TMask, TFeature )
{
  
  std::string localCost = vm["localcost"].as<std::string>();
  std::string initMask = vm["initmask"].as<std::string>();
  std::string outMask = vm["outmask"].as<std::string>();
  std::string initValues = vm["initvalues"].as<std::string>();
  //Output volumes
  std::string output = vm["output"].as<std::string>();
  std::string chosenStcl = vm["chosenstcl"].as<std::string>();
  std::string fmLabels = vm["fmlabels"].as<std::string>();
  //MSFM parameters
  std::string stencilset = vm["stencilset"].as<std::string>();
  std::string msScheme = vm["msscheme"].as<std::string>();
  int fgMask1 = vm["fgmask1"].as<int>();
  int fgMask2 = vm["fgmask2"].as<int>();
  bool dbgMode = vm["debug"].as<bool>();
  bool use26neighs = vm["use26neighs"].as<bool>();
   

//Type definitions
  typedef TFeature	FeaturePixelType;
  typedef TMask		MaskPixelType;
  typedef double	RealPixelType;
  const unsigned int Dimension = 3;
  
  typedef itk::Image<MaskPixelType,Dimension>		MaskImageType;
  typedef itk::Image<FeaturePixelType,Dimension>	FeatureImageType;
  typedef itk::Image<RealPixelType,Dimension>		RealImageType;
  
  typedef itk::ImageFileReader<MaskImageType>		MaskReaderType;
  typedef itk::ImageFileReader<FeatureImageType>	FeatureReaderType;
  typedef itk::ImageFileReader<RealImageType>		RealReaderType;
  
  
  typedef itk::MultiStencilEuclideanFastMarchingImageFilterBase<
		  MaskImageType,FeatureImageType,RealImageType>		FastMarchingBaseType;
  typedef itk::MultiStencilEuclideanFastMarchingCenteredImageFilter<
		    MaskImageType,FeatureImageType,RealImageType>	FMCenteredFilterType;
  typedef itk::MultiStencilEuclideanFastMarching1stOrderImageFilter<
		    MaskImageType,FeatureImageType,RealImageType>	FM1stOrderFilterType;
  typedef itk::MultiStencilEuclideanFastMarching2ndOrderImageFilter<
		    MaskImageType,FeatureImageType,RealImageType>	FM2ndOrderFilterType;
  typedef typename itk::MultiStencilEuclideanFastMarchingImageFilterBase<
	    MaskImageType,FeatureImageType,RealImageType>::Pointer	FastMarchingBasePointer;
  
  typedef typename FastMarchingBaseType::LabelImageType		LabelImageType;
  
  typedef itk::ImageFileWriter<RealImageType>			RealWriterType;
  typedef itk::ImageFileWriter<LabelImageType>			LabelWriterType;
      
//   typedef typename ThicknessFilterType::UseStencilsEnumType	UseStencilsEnumType;
  
  //Read local cost volume
  typename FeatureReaderType::Pointer costReader = FeatureReaderType::New();
  costReader->SetFileName(localCost.c_str());
  
  try { costReader->Update(); std::cout << localCost.c_str() << " opened" << std::endl; } 
  catch ( itk::ExceptionObject & e ) {
    std::cerr << "File " << localCost.c_str() << " could not be read. " << std::endl;
    std::cerr << "Exception caught!" << std::endl << e.GetDescription() << std::endl << e.GetLocation() << std::endl;
    return EXIT_FAILURE;
  }
  
  //Read initial mask
  typename MaskReaderType::Pointer initMaskReader = MaskReaderType::New();
  initMaskReader->SetFileName(initMask.c_str());
  
  try { initMaskReader->Update(); std::cout << initMask.c_str() << " opened" << std::endl; } 
  catch ( itk::ExceptionObject & e ) {
    std::cerr << "File " << initMask.c_str() << " could not be read. " << std::endl;
    std::cerr << "Exception caught!" << std::endl << e.GetDescription() << std::endl << e.GetLocation() << std::endl;
    return EXIT_FAILURE;
  }

  
  //------------------------------- Input Checks -------------------------------
  
  //----------------------- Compute outputs --------------------------------
  
  //SETUP FAST MARCHING FILTER
  FastMarchingBasePointer msfmFilter;
  
  if (msScheme=="MSc")
    msfmFilter = FMCenteredFilterType::New();
  else if (msScheme=="MS1")
    msfmFilter = FM1stOrderFilterType::New();
  else if (msScheme=="MS2")
    msfmFilter = FM2ndOrderFilterType::New();
  else
  {
    std::cerr << "Unknown multi-stencil scheme: " << msScheme << std::endl;
    return EXIT_FAILURE;
  }
    
  
  msfmFilter->SetCostInput( costReader->GetOutput() );
  msfmFilter->SetInitialMaskInput( initMaskReader->GetOutput() );
  msfmFilter->SetMask1Value( (MaskPixelType)fgMask1 );
  msfmFilter->SetMask2Value( (MaskPixelType)fgMask2 );
  msfmFilter->SetUseFullNhood( use26neighs );
  
  if (dbgMode)
    msfmFilter->DebugOn();
  
  
  //Pasar el conjunto de stencils elegido al objeto. Asumimos que la variable stencilset es de tipo std::string
  if (stencilset=="OrthoRSquare1")
    msfmFilter->SetUseStencils( FastMarchingBaseType::OrthoRSquare1 );
  else if (stencilset=="OrthoRSquare2")
    msfmFilter->SetUseStencils( FastMarchingBaseType::OrthoRSquare2 );
  else if (stencilset=="OrthoRSquare5")
    msfmFilter->SetUseStencils( FastMarchingBaseType::OrthoRSquare5 );
  else if (stencilset=="OrthoRSquare6")
    msfmFilter->SetUseStencils( FastMarchingBaseType::OrthoRSquare6 );
  else if (stencilset=="OrthoRSquare9")
    msfmFilter->SetUseStencils( FastMarchingBaseType::OrthoRSquare9 );
  else if (stencilset=="OrthoRSquare10")
    msfmFilter->SetUseStencils( FastMarchingBaseType::OrthoRSquare10 );
  else if (stencilset=="Triangles")
    msfmFilter->SetUseStencils( FastMarchingBaseType::Triangles );
  else if (stencilset=="Hassouna6Stcl")
    msfmFilter->SetUseStencils( FastMarchingBaseType::Hassouna6Stcl );
  else if (stencilset=="HassounaExtended")
    msfmFilter->SetUseStencils( FastMarchingBaseType::HassounaExtended );
  else
  {
    std::cerr << "Unknown stencil set: " << stencilset << " for " << Dimension << "-D images." << std::endl;
    return EXIT_FAILURE;
  }
  
  //------------- Optional input images ----------------------------------------

  if (!outMask.empty() && (outMask.compare("None")!=0) )
  {
    typename MaskReaderType::Pointer outMaskReader = MaskReaderType::New();
    outMaskReader->SetFileName(outMask.c_str());
    
    try { outMaskReader->Update(); std::cout << outMask.c_str() << " opened" << std::endl; } 
    catch ( itk::ExceptionObject & e ) {
      std::cerr << "File " << outMask.c_str() << " could not be read. " << std::endl;
      std::cerr << "Exception caught!" << std::endl << e.GetDescription() << std::endl << e.GetLocation() << std::endl;
      return EXIT_FAILURE;
    }
    
    msfmFilter->SetFinalMaskInput( outMaskReader->GetOutput() );
    msfmFilter->SetMask2Value( (MaskPixelType)fgMask2 );
  }    
  
  
  if (!initValues.empty() && (initValues.compare("None")!=0) )
  {
    typename RealReaderType::Pointer initValuesReader = RealReaderType::New();
    initValuesReader->SetFileName(initValues.c_str());
    
    try { initValuesReader->Update(); std::cout << initValues.c_str() << " opened" << std::endl; } 
    catch ( itk::ExceptionObject & e ) {
      std::cerr << "File " << initValues.c_str() << " could not be read. " << std::endl;
      std::cerr << "Exception caught!" << std::endl << e.GetDescription() << std::endl << e.GetLocation() << std::endl;
      return EXIT_FAILURE;
    }
    
    msfmFilter->SetInitialValuesInput( initValuesReader->GetOutput() );
  }
  
  //-------------------- Write the maps into files -----------------------------
//   try { 
    typename RealWriterType::Pointer writer = RealWriterType::New();
    writer->SetInput( msfmFilter->GetOutput() );
    writer->SetFileName( output.c_str() );
    writer->Update(); 
//   } 
//   catch ( itk::ExceptionObject & e ) {
//     std::cerr << "File " << output.c_str() << " could not be written. " << std::endl;
//     std::cerr << "Exception caught!" << std::endl 
// 	      << e.GetDescription() << std::endl << e.GetLocation() << std::endl;
//     return EXIT_FAILURE;
//   }
    
  if (!chosenStcl.empty() && (chosenStcl.compare("None")!=0) )
  {
    try { 
      typename LabelWriterType::Pointer writerField = LabelWriterType::New();
      writerField->SetInput( msfmFilter->GetChosenStencilOutput() );
      writerField->SetFileName( chosenStcl.c_str() );
      writerField->Update();
    } 
    catch ( itk::ExceptionObject & e ) {
      std::cerr << "File " << chosenStcl.c_str() << " could not be written. " << std::endl;
      std::cerr << "Exception caught!" << std::endl 
      << e.GetDescription() << std::endl << e.GetLocation() << std::endl;
      return EXIT_FAILURE;
    }
  }    

  if ( !fmLabels.empty() && (fmLabels.compare("None")!=0) )
  {
    try { 
      typename LabelWriterType::Pointer writerLabels = LabelWriterType::New();
      writerLabels->SetInput( msfmFilter->GetLabelImage() );
      writerLabels->SetFileName( fmLabels.c_str() );
      writerLabels->Update();
    } 
    catch ( itk::ExceptionObject & e ) {
      std::cerr << "File " << fmLabels.c_str() << " could not be written. " << std::endl;
      std::cerr << "Exception caught!" << std::endl 
      << e.GetDescription() << std::endl << e.GetLocation() << std::endl;
      return EXIT_FAILURE;
    }
  }

  return EXIT_SUCCESS;
}


template<class TMask> 
int DoItAux( po::variables_map &vm, TMask dummy)
{

  //Get the local cost argument
  std::string localCost = vm["localcost"].as<std::string>();
  
  itk::ImageIOBase::IOPixelType pixelType;
  itk::ImageIOBase::IOComponentType componentType;

  try
  {
    itk::GetImageType (localCost, pixelType, componentType);

    // This filter handles all types for the feature pixel.
    
    switch (componentType)
      {
      case itk::ImageIOBase::UCHAR:
        return DoIt( vm, dummy, static_cast<unsigned char>(0));
        break;
      case itk::ImageIOBase::CHAR:
        return DoIt( vm, dummy, static_cast<char>(0));
        break;
      case itk::ImageIOBase::USHORT:
        return DoIt( vm, dummy, static_cast<unsigned short>(0));
        break;
      case itk::ImageIOBase::SHORT:
        return DoIt( vm, dummy, static_cast<short>(0));
        break;
      case itk::ImageIOBase::UINT:
        return DoIt( vm, dummy, static_cast<unsigned int>(0));
        break;
      case itk::ImageIOBase::INT:
        return DoIt( vm, dummy, static_cast<int>(0));
        break;
      case itk::ImageIOBase::ULONG:
        return DoIt( vm, dummy, static_cast<unsigned long>(0));
        break;
      case itk::ImageIOBase::LONG:
        return DoIt( vm, dummy, static_cast<long>(0));
        break;
      case itk::ImageIOBase::FLOAT:
        return DoIt( vm, dummy, static_cast<float>(0));
        break;
      case itk::ImageIOBase::DOUBLE:
        return DoIt( vm, dummy, static_cast<double>(0));
        break;
      case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
      default:
        std::cout << "unknown feature component type" << std::endl;
        break;
      }
    }
  catch( itk::ExceptionObject &excep)
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
  
}



int main( int argc, char * argv[] )
{
  
  po::options_description desc("Allowed options");

  try
  {
    // Manage the input arguments with boost::programoptions
    desc.add_options()
	("help", "produce help message")
	("localcost", po::value<std::string>()->required(), "Local Cost")
	("initmask",  po::value<std::string>()->required(), "Initial Mask Volume")
	("outmask",   po::value<std::string>()->default_value(""), "Output Mask Volume")
	("output",    po::value<std::string>()->required(), "Output")
	("chosenstcl",po::value<std::string>()->default_value(""), "Chosen Stencils")
	("fmlabels",  po::value<std::string>()->default_value(""), "FM final labels")
	("initvalues",po::value<std::string>()->default_value(""), "Initial Values Volume")
	("msscheme",  po::value<std::string>()->default_value("MSc"), "Multi-Stencil scheme")
	("stencilset",po::value<std::string>()->default_value("OrthoRSquare2"), "Stencil Set")
	("fgmask1",   po::value<int>()->default_value(1), "Fg Value Inital Mask")
	("fgmask2",   po::value<int>()->default_value(1), "Fg Value Output Mask")
	("use26neighs", po::value<bool>()->default_value(false), "Use 26 neighbors to expand the front")
	("debug",     po::value<bool>()->default_value(false), "Show debugging output")
    ;

    po::variables_map vm;        
    po::store(po::parse_command_line(argc, argv, desc), vm);
    
    if (vm.count("help"))
    {
      std::cout << "Basic Command Line Parameter App" << std::endl 
                << desc << std::endl; 
      return EXIT_SUCCESS;
    }
    
    po::notify(vm);    
  
    
    itk::ImageIOBase::IOPixelType pixelType;
    itk::ImageIOBase::IOComponentType componentType;
  
    std::string initMask = vm["initmask"].as<std::string>();
    itk::GetImageType (initMask, pixelType, componentType);


    // This filter handles all types
    
    switch (componentType)
      {
      case itk::ImageIOBase::UCHAR:
        return DoItAux( vm, static_cast<unsigned char>(0));
        break;
      case itk::ImageIOBase::CHAR:
        return DoItAux( vm, static_cast<char>(0));
        break;
      case itk::ImageIOBase::USHORT:
        return DoItAux( vm, static_cast<unsigned short>(0));
        break;
      case itk::ImageIOBase::SHORT:
        return DoItAux( vm, static_cast<short>(0));
        break;
      case itk::ImageIOBase::UINT:
        return DoItAux( vm, static_cast<unsigned int>(0));
        break;
      case itk::ImageIOBase::INT:
        return DoItAux( vm, static_cast<int>(0));
        break;
      case itk::ImageIOBase::ULONG:
        return DoItAux( vm, static_cast<unsigned long>(0));
        break;
      case itk::ImageIOBase::LONG:
        return DoItAux( vm, static_cast<long>(0));
        break;
      case itk::ImageIOBase::FLOAT:
        return DoItAux( vm, static_cast<float>(0));
        break;
      case itk::ImageIOBase::DOUBLE:
        return DoItAux( vm, static_cast<double>(0));
        break;
      case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
      default:
        std::cout << "unknown component type" << std::endl;
        break;
      }
    }
  catch( itk::ExceptionObject &excep)
  {
    std::cerr << argv[0] << ": exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
  }
  catch( po::error &e )
  {
    std::cerr << argv[0] << ": exception caught !" << std::endl << e.what() << std::endl;
    std::cerr << desc << std::endl;
    return EXIT_FAILURE;
  }
  
  return EXIT_SUCCESS;
}
