/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2023 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitpit.
 *
 *  bitpit is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  bitpit is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with bitpit. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/

/**
 * \example POD_example_00004.cpp
 * 
 * \brief POD basis computation using voloctree.
 *
 * This example computes the POD basis starting from a database of simulations
 * defined on the same mesh and reconstructs the first mode as a PODField object
 * using the buildFieldWithCoeff function.
 * Following the computation, both the POD Field with the mode and POD Mode object
 * corresponding to the first mode are written on file in the example folder.
 *
 * <b>To run</b>: ./POD_example_00004 \n
 */ 

#include <array>
#if BITPIT_ENABLE_MPI
#include <mpi.h>
#endif

#include "pod.hpp"

using namespace bitpit;

/*
 * Print a matrix.
 *
 * \param[in] mat, matrix of size M rows times N columns.
 */
void printMat (std::vector < std::vector<double>> mat)
{
    std::cout << "mat = " << std::endl;
    size_t M = mat.size();
    size_t N = mat[0].size();
    for (size_t i=0; i<M; i++) {
        for (size_t j=0; j<N; j++) {
            if (j == 0) {
                std::cout << "[ "<< std::setprecision(4) << mat[i][j] ;
            }
            else if (j==(N-1)) {
                std::cout << " , "  << std::setprecision(4) << mat[i][j] << " ]" << std::endl;
            }
            else {
                std::cout << " , " << std::setprecision(4) << mat[i][j] ;
            }
            if (N==1) {
                std::cout  << " ]" << std::setprecision(4) << std::endl;
            }
        }
    }
}

/**
 * Run the example.
 */ 
void run()
{
    /**<Create POD object.*/
    POD pod;

    /**<Add snapshots to database.*/
    pod.addSnapshot("./data", "test_set2.0");
    pod.addSnapshot("./data", "test_set2.1");
    pod.addSnapshot("./data", "test_set2.2");
    pod.addSnapshot("./data", "test_set2.3");
    pod.addSnapshot("./data", "test_set2.4");
    pod.addSnapshot("./data", "test_set2.5");


    /**<Set POD.*/
    pod.setMeshType(POD::MeshType::VOLOCTREE);
    pod.setStaticMesh(true);
    pod.setErrorMode(POD::ErrorMode::NONE);
    pod.setWriteMode(POD::WriteMode::DEBUG);
    pod.setMemoryMode(POD::MemoryMode::MEMORY_NORMAL);
    pod.setEnergyLevel(100.00);
    pod.setUseMean(false);
    pod.setDirectory("pod");
    pod.setName("pod.test.solver");
    pod.setModeCount(6);

    /**<Compute the POD modes.*/
    pod.run();

    std::cout << "the number of modes is = " << pod.getModeCount() << std::endl;
    std::size_t nModes = pod.getModeCount();
    std::size_t nFields = pod.getScalarNames().size()+pod.getVectorNames().size();


    /*
     * Matrici con i prodotti scalari dei modi
     */
    // lista id attivi della POD
    const std::unordered_set<long> & listActIds = pod.getListActiveIDs();
    // modi
    const std::vector<pod::PODMode> &podModes = pod.getModes();
    //mesh
    const VolumeKernel* podMesh = pod.getMesh();

    std::vector<std::vector<double>> prodMatP;
    prodMatP.resize(nModes, std::vector<double>(nModes, 0.0));
    std::vector<std::vector<double>> prodMatT;
    prodMatT.resize(nModes, std::vector<double>(nModes, 0.0));
    std::vector<std::vector<double>> prodMatU;
    prodMatU.resize(nModes, std::vector<double>(nModes, 0.0));

    for (int i=0; i<nModes; i++) {
        for (int j=0; j<nModes; j++) {
            for (long id : listActIds) {
                double modePi = podModes[i].scalar->at(id, 0);
                double modePj = podModes[j].scalar->at(id, 0);
                prodMatP[i][j] += modePi*modePj*podMesh->evalCellVolume(id);

                double modeTi = podModes[i].scalar->at(id, 1);
                double modeTj = podModes[j].scalar->at(id, 1);
                prodMatT[i][j] += modeTi*modeTj*podMesh->evalCellVolume(id);

                std::array<double,3> modeUi = podModes[i].vector->at(id, 0);
                std::array<double,3> modeUj = podModes[j].vector->at(id, 0);
                prodMatU[i][j] += dotProduct(modeUi,modeUj)*podMesh->evalCellVolume(id);
            }
        }
    }
    std::cout<< "prodotti scalari dei modi di pressione" << std::endl;
    printMat(prodMatP);
    std::cout<< "prodotti scalari dei modi di temperatura" << std::endl;
    printMat(prodMatT);
    std::cout<< "prodotti scalari dei modi di velocità" << std::endl;
    printMat(prodMatU);
}

/**
 * Main program.
 */

int main(int argc, char *argv[]) 
{
#if BITPIT_ENABLE_MPI
    MPI_Init(&argc,&argv);
#endif    

    /** Run the example */
    try {
        run();
    } catch (const std::exception &exception) {
        log::cout() << exception.what();
        exit(1);
    }

#if BITPIT_ENABLE_MPI
    MPI_Finalize();
#endif

}
