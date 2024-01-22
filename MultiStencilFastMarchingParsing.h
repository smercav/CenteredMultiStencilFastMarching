
/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $HeadURL: http://svn.slicer.org/Slicer3/branches/Slicer-3-2/Applications/CLI/ProcesadoMRICardiaca/EntradaSalida/CardiacTrasform.cxx $
  Language:  C++
  Date:      $Date: 2006-12-18 15:29:35 +0100 (lun 18 de dic de 2006) $
  Version:   $Revision: 1855 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.
 * 
 * Author: Lucilio Cordero Grande
 * Date: April 15th 2013
 * e-mail: lcorgra@lpi.tel.uva.es

=========================================================================*/

#include <iostream>
#include <iomanip>

#include <typeinfo>
#include <metaCommand.h>

#include "itkCommandLineArgumentParser.h"

struct EstructuraParsing
{
   //Input volumes
   std::string localCost;
   std::string initMask;
   std::string outMask;
   //Output volumes
   std::string output;
   std::string chosenStcl;
   //MSFM parameters
   std::string stencilset;
   int fgMask1;
   int fgMask2;
};

EstructuraParsing parsea(itk::SmartPointer<itk::CommandLineArgumentParser>);

