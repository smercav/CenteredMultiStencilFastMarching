#ifndef __itkPluginUtilities_h
#define __itkPluginUtilities_h

// STD includes
#include <vector>
#include <string>
#include <cstring>

// ITK includes
#include <itkContinuousIndex.h>
#include <itkImage.h>
#include <itkImageFileReader.h>

/**
 * This file is a modification of 3DSlicer itkPluginUtilities.h, 
 * trimmed to provide itk::GetImageType to this project CLIs.
 */

namespace itk
{
  //-----------------------------------------------------------------------------
  /// Get the PixelType and ComponentType from fileName
  void GetImageType (std::string fileName,
                     ImageIOBase::IOPixelType &pixelType,
                     ImageIOBase::IOComponentType &componentType)
    {
      typedef itk::Image<unsigned char, 3> ImageType;
      itk::ImageFileReader<ImageType>::Pointer imageReader =
        itk::ImageFileReader<ImageType>::New();
      imageReader->SetFileName(fileName.c_str());
      imageReader->UpdateOutputInformation();

      pixelType = imageReader->GetImageIO()->GetPixelType();
      componentType = imageReader->GetImageIO()->GetComponentType();
    }
    
  const unsigned int GetImageDimension(std::string fileName)
  {
  
    itk::ImageIOBase::Pointer imageIO =
	  itk::ImageIOFactory::CreateImageIO(
	      fileName.c_str(), itk::ImageIOFactory::ReadMode);
    if( !imageIO )
    {
      std::cerr << "Could not CreateImageIO for: " << fileName << std::endl;
      return EXIT_FAILURE;
    }
    
    imageIO->SetFileName(fileName.c_str());
    imageIO->ReadImageInformation();
    
    return imageIO->GetNumberOfDimensions();
    
  }
//   //-----------------------------------------------------------------------------
//   /// Get the PixelTypes and ComponentTypes from fileNames
//   void GetImageTypes (std::vector<std::string> fileNames,
//                       std::vector<ImageIOBase::IOPixelType> &pixelTypes,
//                       std::vector<ImageIOBase::IOComponentType> &componentTypes)
//     {
//     pixelTypes.clear();
//     componentTypes.clear();
// 
//     // For each file, find the pixel and component type
//     for (std::vector<std::string>::size_type i = 0; i < fileNames.size(); i++)
//       {
//       ImageIOBase::IOPixelType pixelType;
//       ImageIOBase::IOComponentType componentType;
// 
//       GetImageType (fileNames[i],
//                     pixelType,
//                     componentType);
//       pixelTypes.push_back(pixelType);
//       componentTypes.push_back(componentType);
//       }
//     }


} // end namespace itk

#endif
