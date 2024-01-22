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
 * Based on: PerfusionQuantificationParsing.cxx

=========================================================================*/

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include <iostream>
#include <iomanip>

#include <typeinfo>
#include <metaCommand.h>

#include "MultiStencilFastMarchingParsing.h"

EstructuraParsing parsea(itk::SmartPointer<itk::CommandLineArgumentParser> parser)
{
   EstructuraParsing datos;
   
   //Input volumes
   datos.localCost = "";
   parser->GetCommandLineArgument("-localcost",datos.localCost);//Input volume
   
   datos.initMask = "";
   parser->GetCommandLineArgument("-initmask",datos.initMask);//Input volume
   
   datos.outMask = "";
   parser->GetCommandLineArgument("-outmask",datos.outMask);//Input volume
   
   //Output volumes
   datos.output = "";
   parser->GetCommandLineArgument("-output",datos.output);//Input volume
   
   datos.chosenStcl = "";
   parser->GetCommandLineArgument("-localcost",datos.chosenStcl);//Input volume
   
   //MSFM parameters
   datos.stencilset = "OrthoRSquare2";
   parser->GetCommandLineArgument("-stencilset", datos.stencilset);
   
   datos.fgMask1=1;
   parser->GetCommandLineArgument("-fgmask1",datos.fgMask1);
   
   datos.fgMask2=1;
   parser->GetCommandLineArgument("-fgmask2",datos.fgMask2);
   
   return datos;
}
