/* libFsGrid
This file is part of Analysator.
Copyright 2021 University of Helsinki

For details of usage, see the COPYING file and read the "Rules of the Road"
at http: //www.physics.helsinki.fi/vlasiator/

This program is free software; you can redistribute it and / or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
 51 Franklin Street, Fifth Floor, Boston, MA 02110 - 1301 USA.
*/

#include <iostream>
#include <stdlib.h>
#include <vector>
#include <vlsv_reader.h>
#include <vlsv_common.h>
#include <cmath>
#include <cstring>

using namespace std;
typedef double Real;
bool readFsGridVariable(std::string, std::string, int c, std::vector<Real> &rawData, std::vector<int> &shape);

#ifdef __cplusplus
extern "C"
#endif

//Wrapper entry point to interface with python
bool read(char *fname,char *varName ,int component,Real *vec,int* ret_shape){

   std::vector<Real> rawData;
   std::vector<int> shape;
   shape.resize(3);
   if (!readFsGridVariable(fname, varName, component, rawData, shape)){
      return false;
   };
   std::memcpy(vec, rawData.data(), rawData.size() * sizeof(Real));
   std::memcpy(ret_shape, shape.data(), shape.size() * sizeof(int));
   return true;
}


// ! Helper function: calculate position of the local coordinate space for the given dimension
// param globalCells Number of cells in the global Simulation, in this dimension
// param ntasks Total number of tasks in this dimension
// param my_n This task's position in this dimension
// return Cell number at which this task's domains cells start (actual cells, not counting ghost cells)
int32_t calcLocalStart(int32_t globalCells, int ntasks, int my_n){
   int n_per_task = globalCells / ntasks;
   int remainder = globalCells % ntasks;

   if (my_n < remainder)
   {
      return my_n * (n_per_task + 1);
   }
   else
   {
      return my_n * n_per_task + remainder;
   }
}

//! Helper function: calculate size of the local coordinate space for the given dimension
// \param globalCells Number of cells in the global Simulation, in this dimension
// \param ntasks Total number of tasks in this dimension
// \param my_n This task's position in this dimension
// \return Nmuber of cells for this task's local domain (actual cells, not counting ghost cells)
int32_t calcLocalSize(int32_t globalCells, int ntasks, int my_n){
   int n_per_task = globalCells / ntasks;
   int remainder = globalCells % ntasks;
   if (my_n < remainder)
   {
      return n_per_task + 1;
   }
   else
   {
      return n_per_task;
   }
}

//! Helper function to optimize decomposition of this grid over the given number of tasks
void computeDomainDecomposition(const std::array<int, 3> &GlobalSize, int nProcs, std::array<int, 3> &processDomainDecomposition){
   std::array<double, 3> systemDim;
   std::array<double, 3> processBox;
   double optimValue = std::numeric_limits<double>::max();
   for (int i = 0; i < 3; i++)
   {
      systemDim[i] = (double)GlobalSize[i];
   }
   processDomainDecomposition = {1, 1, 1};
   for (int i = 1; i <= std::min(nProcs, GlobalSize[0]); i++)
   {
      processBox[0] = std::max(systemDim[0] / i, 1.0);
      for (int j = 1; j <= std::min(nProcs, GlobalSize[1]); j++)
      {
         if (i * j > nProcs)
            break;
         processBox[1] = std::max(systemDim[1] / j, 1.0);
         for (int k = 1; k <= std::min(nProcs, GlobalSize[2]); k++)
         {
            if (i * j * k > nProcs)
               break;
            processBox[2] = std::max(systemDim[2] / k, 1.0);
            double value =
               10 * processBox[0] * processBox[1] * processBox[2] +
               (i > 1 ? processBox[1] * processBox[2] : 0) +
               (j > 1 ? processBox[0] * processBox[2] : 0) +
               (k > 1 ? processBox[0] * processBox[1] : 0);

            if (value < optimValue)
            {
               optimValue = value;
               processDomainDecomposition[0] = i;
               processDomainDecomposition[1] = j;
               processDomainDecomposition[2] = k;
            }
         }
      }
   }

   if (optimValue == std::numeric_limits<double>::max() ||
      processDomainDecomposition[0] * processDomainDecomposition[1] * processDomainDecomposition[2] != nProcs)
   {
      std::cerr << "FSGrid domain decomposition failed, are you running on a prime number of tasks?" << std::endl;
      throw std::runtime_error("FSGrid computeDomainDecomposition failed");
   }
}

//! Helper function to actually read the FsGrid variable data
bool readFsGridVariable(std::string file,std::string varToExtract,int compToExtract,std::vector<Real> &rawData,std::vector<int> &shape){

   // Get Spatial Grid's  max refinement Level
   int maxRefLevel = 0;
   list<pair<string, string>> meshAttributesIn;
   meshAttributesIn.push_back(make_pair("name", "SpatialGrid"));
   map<string,string> meshAttributesOut;
   vlsv::Reader vlsvReader;
   vlsvReader.open(file);

   if (vlsvReader.getArrayAttributes("MESH", meshAttributesIn,meshAttributesOut) == false)
   {
      return false;
   }

   bool meshSuccess = true;
   bool variableSuccess = true;

   vlsv::datatype::type meshDataType;
   vlsv::datatype::type variableDataType;
   uint64_t meshArraySize, meshVectorSize, meshDataSize;
   uint64_t variableArraySize, variableVectorSize, variableDataSize;

   list<pair<string, string>> variableAttributes;
   variableAttributes.push_back(make_pair("mesh", "fsgrid"));
   variableAttributes.push_back(make_pair("name", varToExtract));
   //Read in array size, vector size, data type and data size of the array "VARIABLE" in the vlsv file (Needed in reading the array)

   if (vlsvReader.getArrayInfo("VARIABLE", variableAttributes, variableArraySize, variableVectorSize, variableDataType, variableDataSize) == false)
   {
      return false;
   }

   std::map<string, string>::iterator attributesOutIt;
   attributesOutIt = meshAttributesOut.find("max_refinement_level");
   if (attributesOutIt != meshAttributesOut.end())
   {
      maxRefLevel = stoi(attributesOutIt->second);
   }
   int numtasks;
   int xcells,ycells,zcells;
   vlsvReader.readParameter("numWritingRanks",numtasks);
   vlsvReader.readParameter("xcells_ini",xcells);
   vlsvReader.readParameter("ycells_ini",ycells);
   vlsvReader.readParameter("zcells_ini",zcells);
   xcells*=pow(2,maxRefLevel);
   ycells*=pow(2,maxRefLevel);
   zcells*=pow(2,maxRefLevel);
   std::array<int,3> GlobalBox={xcells,ycells,zcells};
   std::array<int,3> thisDomainDecomp;

   shape[0]=xcells;
   shape[1]=ycells;
   shape[2]=zcells;
   //Compute Domain Decomposition Scheme for this vlsv file
   computeDomainDecomposition(GlobalBox,numtasks,thisDomainDecomp);
   std::map<uint, Real> orderedData;
   rawData.resize(xcells * ycells * zcells);
   std::array<int32_t,3> taskSize,taskStart;
   std::array<int32_t,3> taskEnd;
   int readOffset=0;
   int index,my_x,my_y,my_z;

   for (int task=0; task<numtasks; task++){

      my_x=task/thisDomainDecomp[2]/thisDomainDecomp[1];
      my_y=(task/thisDomainDecomp[2])%thisDomainDecomp[1];
      my_z=task%thisDomainDecomp[2];

      taskStart[0] = calcLocalStart(GlobalBox[0], thisDomainDecomp[0] ,my_x);
      taskStart[1] = calcLocalStart(GlobalBox[1], thisDomainDecomp[1] ,my_y);
      taskStart[2] = calcLocalStart(GlobalBox[2], thisDomainDecomp[2] ,my_z);

      taskSize[0] = calcLocalSize(GlobalBox[0], thisDomainDecomp[0] ,my_x);
      taskSize[1] = calcLocalSize(GlobalBox[1], thisDomainDecomp[1] ,my_y);
      taskSize[2] = calcLocalSize(GlobalBox[2], thisDomainDecomp[2] ,my_z);

      taskEnd[0]= taskStart[0]+taskSize[0];
      taskEnd[1]= taskStart[1]+taskSize[1];
      taskEnd[2]= taskStart[2]+taskSize[2];

      int64_t readSize=  taskSize[0] * taskSize[1] * taskSize[2] ;
      //Allocate vector for reading
      std::vector<Real> buffer(readSize*variableVectorSize);

      if ( variableDataSize==sizeof(Real)){
         if (vlsvReader.readArray("VARIABLE", variableAttributes, readOffset, readSize,  (char*)buffer.data()) == false) {
            variableSuccess = false;
            return false;
         }
      }else{
         std::vector<float> tmpbuffer(readSize * variableVectorSize);
         if (vlsvReader.readArray("VARIABLE", variableAttributes, readOffset, readSize, (char *)tmpbuffer.data()) == false){
            variableSuccess = false;
            return false;
         }
         for (unsigned int i = 0; i < readSize * variableVectorSize; i++){
            buffer[i] = tmpbuffer[i];
         }
      }

      uint64_t globalindex,counter=0;;
      for (int z=taskStart[2]; z<taskEnd[2]; z++){
         for (int y=taskStart[1]; y< taskEnd[1]; y++){
            for (int x=taskStart[0]; x<taskEnd[0]; x++){
               globalindex= x + y*xcells + z*xcells*ycells;
               switch (variableDataType){
                  case vlsv::datatype::type::FLOAT:
                     if (variableDataSize == sizeof(float))
                        memcpy(&rawData[globalindex], &buffer[counter+compToExtract],sizeof(double));
                     if (variableDataSize == sizeof(double))
                        memcpy(&rawData[globalindex], &buffer[counter + compToExtract],  sizeof(double));
                     break;
                  case vlsv::datatype::type::UINT:
                     memcpy(&rawData[globalindex], &buffer[counter + compToExtract],  sizeof(uint));
                     break;
                  case vlsv::datatype::type::INT:
                     memcpy(&rawData[globalindex], &buffer[counter + compToExtract],  sizeof(int));
                     break;
                  case vlsv::datatype::type::UNKNOWN:
                     cerr << "ERROR, BAD DATATYPE AT " << __FILE__ << " " << __LINE__ << endl;
                     return false;
                     break;
               }
               counter+=variableVectorSize;
            }
         }
      }
      readOffset+=readSize;
//
   }
   return true;
}



