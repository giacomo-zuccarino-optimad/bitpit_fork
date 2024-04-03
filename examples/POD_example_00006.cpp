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
 * Build a PODField object with difference between two PODField objects.
 *
 * \param[in] field1, minuend PODField.
 * \param[in] field2, subtrahend PODField.
 * \param[in] listActIds, list of active ids of the mesh of field1 and field2.
 * \param[out] errorField, PODField with the difference between field1 and field2.
 */
pod::PODField buildErrorField (pod::PODField &field1, pod::PODField &field2, const std::unordered_set<long> listActIds)
{
    pod::PODField errorField;
    errorField.scalar = std::unique_ptr<pod::ScalarStorage>(new pod::ScalarStorage(2, &field1.mesh->getCells()));
    errorField.vector = std::unique_ptr<pod::VectorStorage>(new pod::VectorStorage(1, &field1.mesh->getCells()));
    errorField.scalar->fill(0.0);
    errorField.vector->fill(std::array<double, 3>{{0.0, 0.0, 0.0}});
    errorField.mask = std::unique_ptr<PiercedStorage<bool>>(new PiercedStorage<bool>(1, &field1.mesh->getCells()));
    errorField.mask->fill(0.0);
    errorField.setMesh(field1.mesh);
    for (long id : listActIds) {
        for (std::size_t isf = 0; isf < 2; ++isf) {
            double field1s = field1.scalar->at(id, isf);
            double field2s = field2.scalar->at(id, isf);
            errorField.scalar->at(id, isf) += field1s-field2s;
        }
        std::array<double,3> field1sv = field1.vector->at(id,0);
        std::array<double,3> field2sv = field2.vector->at(id,0);
        errorField.vector->at(id, 0) += field1sv-field2sv;
        errorField.mask->at(id) = 1;
    }
    return errorField;
}

/*
 * Build a PODField object with the relative difference between two PODField objects.
 *
 * \param[in] field1, minuend PODField.
 * \param[in] field2, subtrahend PODField.
 * \param[in] pod, POD object over which the two fields have been defined.
 * \param[in] norm_type, L2 or Linf
 * \param[out] errorField, PODField with the difference between field1 and field2.
 */
pod::PODField buildRelErrorField (pod::PODField &field1, pod::PODField &field2, POD &pod, std::string norm_type )
{
    pod::PODField errorField;
    errorField.scalar = std::unique_ptr<pod::ScalarStorage>(new pod::ScalarStorage(2, &field1.mesh->getCells()));
    errorField.vector = std::unique_ptr<pod::VectorStorage>(new pod::VectorStorage(1, &field1.mesh->getCells()));
    errorField.scalar->fill(0.0);
    errorField.vector->fill(std::array<double, 3>{{0.0, 0.0, 0.0}});
    errorField.mask = std::unique_ptr<PiercedStorage<bool>>(new PiercedStorage<bool>(1, &field1.mesh->getCells()));
    errorField.mask->fill(0.0);
    errorField.setMesh(field1.mesh);
    std::vector<double> vec;
    if (norm_type == "L2") {
        vec = pod.fieldsl2norm(field1);
    }
    else if (norm_type == "Linf") {
        vec = pod.fieldsMax(field1);
    }
    const std::unordered_set<long> & listActIds = pod.getListActiveIDs();
    for (long id : listActIds) {
        for (std::size_t isf = 0; isf < 2; ++isf) {
            double field1s = field1.scalar->at(id, isf);
            double field2s = field2.scalar->at(id, isf);
            errorField.scalar->at(id, isf) += (field1s-field2s)/vec[isf];
        }
        std::array<double,3> field1sv = field1.vector->at(id,0);
        std::array<double,3> field2sv = field2.vector->at(id,0);
        errorField.vector->at(id, 0) += (field1sv-field2sv)/vec[2];
        errorField.mask->at(id) = 1;
    }
    return errorField;
}


/*
 * Print the L2 norm of each field of a PODField object.
 *
 * \param[in] field, PODField object.
 * \param[in] pod, POD object defined on the same mesh.
 * \param[in] field_name, string of the field name to display.
 */
void printL2norm (pod::PODField &field, POD &pod, std::string field_name)
{
    std::vector<double> vecL2 = pod.fieldsl2norm(field);
    std::vector<std::string> scalarNames = pod.getScalarNames();
    std::vector<std::array<std::string,3>> vectorNames= pod.getVectorNames();
    int N = scalarNames.size();
    for (int i=0; i<N; i++) {
        //std::cout << "L2 norm of " << field_name << " " << scalarNames[i] << " is "<< vecL2[i] << std::endl;
    }
    int M = vectorNames.size();
    for (int i=N; i<N+M; i++) {
        //std::cout << "L2 norm of " << field_name << " " << vectorNames[i-N][0].substr(0,vectorNames[i-N][0].size()-2) << " is "<< vecL2[i] << std::endl;
    }
    double norm = std::sqrt(dotProduct(vecL2,vecL2));
    std::cout << "total L2 norm of " << field_name << " is "<< norm << std::endl;
}

/*
 * Print the L infinity norm of each field of a PODField object.
 *
 * \param[in] field, PODField object.
 * \param[in] pod, POD object defined on the same mesh.
 * \param[in] field_name, string of the field name to display.
 */
void printLinfnorm (pod::PODField &field, POD &pod, std::string field_name)
{
    std::vector<double> vecLinf = pod.fieldsMax(field);
    std::vector<std::string> scalarNames = pod.getScalarNames();
    std::vector<std::array<std::string,3>> vectorNames= pod.getVectorNames();
    int N = scalarNames.size();
    for (int i=0; i<N; i++) {
        std::cout << "L infinity norm of " << field_name << " " << scalarNames[i] << " is "<< vecLinf[i] << std::endl;
    }
    int M = vectorNames.size();
    for (int i=N; i<N+M; i++) {
        std::cout << "L infinity norm of " << field_name << " " << vectorNames[i-N][0].substr(0,vectorNames[i-N][0].size()-2) << " is "<< vecLinf[i] << std::endl;
    }
}

void printMatFinal (std::vector < std::vector<double>> mat)
{
    std::cout << "mat = " << std::endl;
    size_t M = mat.size();
    size_t N = mat[0].size();
    for (size_t i=0; i<M; i++) {
        for (size_t j=0; j<N; j++) {
            if (j == 0) {
                std::cout << "[ "<< std::setprecision(12) << mat[i][j] ;
            }
            else if (j==(N-1)) {
                std::cout << " , "  << std::setprecision(12) << mat[i][j] << " ]" << std::endl;
            }
            else {
                std::cout << " , " << std::setprecision(12) << mat[i][j] ;
            }
            if (N==1) {
                std::cout  << " ]" << std::endl;
            }
        }
    }
}

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
                std::cout << "[ "<< mat[i][j] ;
            }
            else if (j==(N-1)) {
                std::cout << " , "  << mat[i][j] << " ]" << std::endl;
            }
            else {
                std::cout << " , " << mat[i][j] ;
            }
            if (N==1) {
                std::cout  << " ]" << std::endl;
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
    /*
    for (int i=0; i<2; i++) {
        pod.addSnapshot("./data", "test_set2."+std::to_string(i));
    }
    */
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

    /**<Compute the POD basis.*/
    pod.compute();

    /* Remark: the reconstruction of the first mode is equal to the first mode
     * only if the mean is not used in the POD algorithm
     * if the mean is used, the mean has to be subtracted to the reconstruction
     * to get the first mode */

    /**<Reconstruc the first mode. */
    std::cout << "the number of modes is = " << pod.getModeCount() << std::endl;
    std::size_t N_modes = pod.getModeCount();
    std::size_t N_fields = pod.getScalarNames().size()+pod.getVectorNames().size();

    /* set up of the coefficient matrix
     * each column contains the coefficient of a specific mode
     * each row contains the coefficient of a specific field, either scalar or vector.*/
    std::vector < std::vector<double>> coeff_mat;
    coeff_mat.clear();
    coeff_mat.resize(N_fields, std::vector<double>(N_modes, 0.0));
    for (std::size_t i = 0; i < N_fields; i++) {
        coeff_mat[i][0] = 1;
    }
    pod::PODField mode_1;
    pod.buildFieldsWithCoeff(coeff_mat, mode_1);
    printL2norm(mode_1, pod, "mode 1");
    std::vector < std::vector<double>> mode1projMat = pod.projectField(mode_1);
    std::cout << "the coefficient matrix of the projection on the first mode is " << std::endl;
    printMat(mode1projMat);

    /* <Write the first mode as VTK file */
    pod.write(mode_1, "mode.1");
    pod.write(0, "mode.1");
    pod.dumpField("mode.1.dump",mode_1);

    pod::PODField mode_1_bwc; // bwc = built with coefficients
    // letto da file perchè la ricostruzione accetta file letti
    pod::SnapshotFile snap ("pod", "mode.1.dump");
    pod.readSnapshot(snap,mode_1_bwc);

    pod::PODField mode_1_recon; // recon = reconstructed with least square
    pod.reconstructFields(mode_1_bwc, mode_1_recon);
    std::vector < std::vector<double>> mode1reconMat = pod.getReconstructionCoeffs();
    printMat(mode1reconMat);

    // estrazione matrice autovalori
    std::vector < std::vector<double>> lambdaMat = pod.getLambda();
    printMatFinal(lambdaMat);

    const std::unordered_set<long> & listActIds = pod.getListActiveIDs();

    //std::vector<double> l2uvec(6,0.0);
    for (int i=0; i<6; i++) {
        std::cout << "confronto per lo snapshot numero == " << i << std::endl;
        const pod::SnapshotFile snap_i ("./data", "test_set2."+std::to_string(i));
        pod::PODField field_i, field_i_recon;
        pod.readSnapshot(snap_i,field_i);
        pod.reconstructFields(field_i, field_i_recon);
        pod::PODField error_i = buildErrorField(field_i, field_i_recon, listActIds);
        pod::PODField error_i_relL2 = buildRelErrorField(field_i, field_i_recon, pod, "L2");
        //printL2norm(error_i,pod,"reconstruction error");
        printL2norm(error_i_relL2,pod,"reconstruction error (relative)");
        //printL2norm(field_i,pod, "field norm");
    }

    for (int i=0; i<6; i++) {
        std::cout << "confronto per lo snapshot numero == " << i << std::endl;
        const pod::SnapshotFile snap_i ("./data", "test_set2."+std::to_string(i));
        pod::PODField field_i, proj_i;
        pod.readSnapshot(snap_i,field_i);
        std::vector<std::vector<double> > projicoeffMat = pod.projectField(field_i);
        pod.buildFieldsWithCoeff(projicoeffMat, proj_i);
        pod::PODField error_ip = buildErrorField(field_i, proj_i, listActIds);
        pod::PODField error_ip_relL2 = buildRelErrorField(field_i, proj_i, pod, "L2");
        //printL2norm(error_ip,pod,"projection error");
        printL2norm(error_ip_relL2,pod,"projection error (relative)");
    }

    /*
     * ricostruzione e confronto con uno snapshot diverso da un generatore
     */

    std::cout << "confronto con angolo di attacco alpha = 1.5 " << std::endl;

    const pod::SnapshotFile snap_15 ("./data", "test_set2.1.5");
    pod::PODField field_15, field_15_recon;
    pod.readSnapshot(snap_15,field_15);
    pod.reconstructFields(field_15, field_15_recon);

    pod::PODField error15 = buildErrorField(field_15, field_15_recon, listActIds);
    pod::PODField error15relL2 = buildRelErrorField(field_15, field_15_recon, pod, "L2");
    //printL2norm(error15,pod,"reconstruction error");
    printL2norm(error15relL2,pod,"reconstruction error (relative)");

    // confronto tra matrici di ricostruzione
    std::vector<std::vector<double> > rec15coeffMat = pod.getReconstructionCoeffs();
    //printMatFinal(rec15coeffMat);
    std::vector<std::vector<double> > proj15coeffMat = pod.projectField(field_15);
    //printMatFinal(proj15coeffMat);
    rec15coeffMat -= proj15coeffMat;
    std::cout << "differenza ricostruzione-proiezione " << std::endl;
    printMatFinal(rec15coeffMat);

    //confronto con la proiezione
    pod::PODField proj_15;
    pod.buildFieldsWithCoeff(proj15coeffMat, proj_15);
    pod::PODField error15p = buildErrorField(field_15, proj_15, listActIds);
    pod::PODField error15prelL2 = buildRelErrorField(field_15, proj_15, pod, "L2");
    //printL2norm(error15p,pod,"projection error");
    printL2norm(error15prelL2,pod,"projection error (relative)");

    std::cout << "confronto con angolo di attacco alpha = 4.5 " << std::endl;

    const pod::SnapshotFile snap_45 ("./data", "test_set2.4.5");
    pod::PODField field_45, field_45_recon;
    pod.readSnapshot(snap_45,field_45);
    pod.reconstructFields(field_45, field_45_recon);

    pod::PODField error45 = buildErrorField(field_45, field_45_recon, listActIds);
    pod::PODField error45relL2 = buildRelErrorField(field_45, field_45_recon, pod, "L2");
    //printL2norm(error45,pod,"reconstruction error");
    printL2norm(error45relL2,pod,"reconstruction error (relative)");

    // confronto tra matrici di ricostruzione
    std::vector<std::vector<double> > rec45coeffMat = pod.getReconstructionCoeffs();
    //printMatFinal(rec45coeffMat);
    std::vector<std::vector<double> > proj45coeffMat = pod.projectField(field_45);
    //printMatFinal(proj45coeffMat);
    rec45coeffMat -= proj45coeffMat;
    std::cout << "differenza ricostruzione-proiezione " << std::endl;
    printMatFinal(rec45coeffMat);

    //confronto con la proiezione
    pod::PODField proj_45;
    pod.buildFieldsWithCoeff(proj45coeffMat, proj_45);
    pod::PODField error45p = buildErrorField(field_45, proj_45, listActIds);
    pod::PODField error45prelL2 = buildRelErrorField(field_45, proj_45, pod, "L2");
    //printL2norm(error45p,pod,"projection error");
    printL2norm(error45prelL2,pod,"projection error (relative)");

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
