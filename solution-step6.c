// Translate this file with
//
// g++ -O3 --std=c++11 assignment-2019.c -o assignment
//
// Run it with
//
// ./assignment
//
// There should be a result.pvd file that you can open with Paraview.
// Sometimes, Paraview requires to select the representation "Point Gaussian"
// to see something meaningful.
//
// (C) 2018-2019 Tobias Weinzierl

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <math.h>
#include <limits>
#include <iomanip>
#include <omp.h>

double t = 0;
double tFinal = 0;
double tPlot = 0;
double tPlotDelta = 0;

int NumberOfBodies = 0;

/**
 * Pointer to pointers. Each pointer in turn points to three coordinates, i.e.
 * each pointer represents one molecule/particle/body.
 */
double** x;

/**
 * Equivalent to x storing the velocities.
 */
double** v;

/**
 * One mass entry per molecule/particle.
 */
double* mass;

/**
 * Global time step size used.
 */
double   timeStepSize = 0.0;

/**
 * Maximum velocity of all particles.
 */
double   maxV;

/**
 * Minimum distance between two elements.
 */
double   minDx;

/**
 * Set up scenario from the command line.
 *
 * This operation is not to be changed in the assignment.
 */
void setUp(int argc, char** argv) {
    NumberOfBodies = (argc - 4) / 7;

    x = new double* [NumberOfBodies];
    v = new double* [NumberOfBodies];
    mass = new double[NumberOfBodies];

    int readArgument = 1;

    tPlotDelta = std::stof(argv[readArgument]); readArgument++;
    tFinal = std::stof(argv[readArgument]); readArgument++;
    timeStepSize = std::stof(argv[readArgument]); readArgument++;

    for (int i = 0; i < NumberOfBodies; i++) {
        x[i] = new double[3];
        v[i] = new double[3];

        x[i][0] = std::stof(argv[readArgument]); readArgument++;
        x[i][1] = std::stof(argv[readArgument]); readArgument++;
        x[i][2] = std::stof(argv[readArgument]); readArgument++;

        v[i][0] = std::stof(argv[readArgument]); readArgument++;
        v[i][1] = std::stof(argv[readArgument]); readArgument++;
        v[i][2] = std::stof(argv[readArgument]); readArgument++;

        mass[i] = std::stof(argv[readArgument]); readArgument++;

        if (mass[i] <= 0.0) {
            std::cerr << "invalid mass for body " << i << std::endl;
            exit(-2);
        }
    }

    std::cout << "created setup with " << NumberOfBodies << " bodies" << std::endl;

    if (tPlotDelta <= 0.0) {
        std::cout << "plotting switched off" << std::endl;
        tPlot = tFinal + 1.0;
    }
    else {
        std::cout << "plot initial setup plus every " << tPlotDelta << " time units" << std::endl;
        tPlot = 0.0;
    }
}


std::ofstream videoFile;


/**
 * This operation is not to be changed in the assignment.
 */
void openParaviewVideoFile() {
    videoFile.open("result.pvd");
    videoFile << "<?xml version=\"1.0\"?>" << std::endl
        << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl
        << "<Collection>";
}





/**
 * This operation is not to be changed in the assignment.
 */
void closeParaviewVideoFile() {
    videoFile << "</Collection>"
        << "</VTKFile>" << std::endl;
}


/**
 * The file format is documented at http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
 *
 * This operation is not to be changed in the assignment.
 */
void printParaviewSnapshot() {
    static int counter = -1;
    counter++;
    std::stringstream filename;
    filename << "result-" << counter << ".vtp";
    std::ofstream out(filename.str().c_str());
    out << "<VTKFile type=\"PolyData\" >" << std::endl
        << "<PolyData>" << std::endl
        << " <Piece NumberOfPoints=\"" << NumberOfBodies << "\">" << std::endl
        << "  <Points>" << std::endl
        << "   <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">";
    //      << "   <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">";

    for (int i = 0; i < NumberOfBodies; i++) {
        out << x[i][0]
            << " "
            << x[i][1]
            << " "
            << x[i][2]
            << " ";
    }

    out << "   </DataArray>" << std::endl
        << "  </Points>" << std::endl
        << " </Piece>" << std::endl
        << "</PolyData>" << std::endl
        << "</VTKFile>" << std::endl;

    videoFile << "<DataSet timestep=\"" << counter << "\" group=\"\" part=\"0\" file=\"" << filename.str() << "\"/>" << std::endl;
}



/**
 * This is the only operation you are allowed to change in the assignment.
 */
void updateBody() {


    maxV = 0.0;
    minDx = std::numeric_limits<double>::max();

    double* force0 = new double[NumberOfBodies];
    double* force1 = new double[NumberOfBodies];
    double* force2 = new double[NumberOfBodies];
    double* forces = new double[NumberOfBodies * 3];
    double* velocities = new double[NumberOfBodies];
    int numBuckets = 10;
    int** buckets = new int* [numBuckets];
    for (int i = 0; i < numBuckets; ++i) {
        buckets[i] = new int[NumberOfBodies + 1];
    }
    double vbucket = 0.0;

    for (int i = 0; i < numBuckets; i++) {
        buckets[i][0] = 0.0;
    }

    // storing the velocities and max velocity of the particles
    // Parallelising this loop does not improve performance as there is very little work done in each thread.
    // The overhead of setting up the threads outweighs the benefits of parallelism and makes the code run slower
    for (int i = 0; i < NumberOfBodies; i++) {
        velocities[i] = sqrt(v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]);
        if (i == 0) {
            maxV = std::sqrt(v[0][0] * v[0][0] + v[0][1] * v[0][1] + v[0][2] * v[0][2]);
        }
        else {
            maxV = std::max(maxV, sqrt(v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]));
        }
    }


    vbucket = maxV / numBuckets;

    // Parallelising this loop does not improve performance as there is very little work done in each thread.
    // The overhead of setting up the threads outweighs the benefits of parallelism and makes the code run slower
    for (int i = 0; i < NumberOfBodies; i++) {
        //printf("Loop: %u Thread: %u\n", i, omp_get_thread_num());
        bool assigned = false;

        for (int j = 0; j < numBuckets; j++) {
            //adding an element to the bucket if its velocity is bounded by the bucket
            if ((j)*vbucket <= velocities[i] && velocities[i] < ((j + 1) * vbucket)) {
                buckets[j][int(buckets[j][0]) + 1] = i;
                buckets[j][0] += 1;//increment the counter of objects in the bucket
                assigned = true;
                break;
            }
        }
        if (assigned == false) { //must be maxV and in last bucket
            buckets[numBuckets - 1][int(buckets[numBuckets - 1][0]) + 1] = i;
            buckets[numBuckets - 1][0] += 1;
        }
    }


    // Create arrays to store collisions within the main loop

    // Collisions will be processed after the main loop has finished
    int collisionTotal = 0;
    int* collision1 = new int[NumberOfBodies];
    int* collision2 = new int[NumberOfBodies];

    // In the last bucket each body will undergo 2^numbuckets microSteps in each updateBody()
    int nummicrosteps = pow(2, numBuckets);


    // The code runs around 3x as fast when compiled with -fopenmp.

    // Run each microstep in parallel on as many threads as possible.
#pragma omp parallel for reduction(min:minDx)
    for (int microstep = 0; microstep < nummicrosteps; microstep++) {

        //loop through the buckets
        for (int i = 0; i < numBuckets; i++) {

            int runEvery = int(pow(2, (numBuckets-1) - i));

            // Only process the bucket if it is due to be processed:
            // i.e. Process:
            //  The last bucket on every microstep
            //  The 2nd last bucket on every 2nd microstep
            //  The 3rd last bucket on every 4th microstep, etc

            if (microstep % runEvery == 0) {

                // Calculate the micro timestep
                float newTimeStepSize = runEvery * timeStepSize / nummicrosteps;

                for (int k = 1; k < buckets[i][0] + 1; k++) { //for all objects in bucket i

                    force0[buckets[i][k]] = 0.0;
                    force1[buckets[i][k]] = 0.0;
                    force2[buckets[i][k]] = 0.0;

                    for (int l = 0; l < NumberOfBodies; l++) {

                        if (buckets[i][k] != l) {
                            const double distance = sqrt(
                                (x[buckets[i][k]][0] - x[l][0]) * (x[buckets[i][k]][0] - x[l][0]) +
                                (x[buckets[i][k]][1] - x[l][1]) * (x[buckets[i][k]][1] - x[l][1]) +
                                (x[buckets[i][k]][2] - x[l][2]) * (x[buckets[i][k]][2] - x[l][2])
                            );

                            // x,y,z forces acting on particle 0
                            force0[buckets[i][k]] += (x[l][0] - x[buckets[i][k]][0]) * mass[l] * mass[buckets[i][k]] / distance / distance / distance;
                            force1[buckets[i][k]] += (x[l][1] - x[buckets[i][k]][1]) * mass[l] * mass[buckets[i][k]] / distance / distance / distance;
                            force2[buckets[i][k]] += (x[l][2] - x[buckets[i][k]][2]) * mass[l] * mass[buckets[i][k]] / distance / distance / distance;

                            if (k < l) {// no need to calculate twice
                                minDx = std::min(minDx, distance);
                                //collision if distance between 2 objects less than sum of radii
                                if (distance <= 0.01) {
                                #pragma opemmp critical
                                    {
                                        // Add these 2 bodies to the collision list
                                        // Do not add to the collision list if this collision has already been noted
                                        bool found = false;
                                        for (int ic = 0; ic < collisionTotal; ic++) {
                                            if (((l == collision1[ic]) && (buckets[i][k] == collision2[ic])) || ((buckets[i][k] == collision1[ic]) && (l == collision2[ic]))) {
                                                found = true;
                                            }
                                        }
                                        if (!found) {
                                            collision1[collisionTotal] = l;
                                            collision2[collisionTotal] = buckets[i][k];
                                            collisionTotal++;
                                        }
                                    }


                                }
                            }
                        }
                    }
                    //storing the forces acting on each object
                    forces[3 * buckets[i][k]] = force0[buckets[i][k]];
                    forces[(3 * buckets[i][k]) + 1] = force1[buckets[i][k]];
                    forces[(3 * buckets[i][k]) + 2] = force2[buckets[i][k]];
                }
                //updating the positions, velocities, maxV and mass depending on if the particle collided
                for (int k = 1; k < buckets[i][0] + 1; k++) {
                        x[buckets[i][k]][0] = x[buckets[i][k]][0] + newTimeStepSize * v[buckets[i][k]][0];
                        x[buckets[i][k]][1] = x[buckets[i][k]][1] + newTimeStepSize * v[buckets[i][k]][1];
                        x[buckets[i][k]][2] = x[buckets[i][k]][2] + newTimeStepSize * v[buckets[i][k]][2];

                        v[buckets[i][k]][0] = v[buckets[i][k]][0] + newTimeStepSize * force0[buckets[i][k]] / mass[buckets[i][k]];
                        v[buckets[i][k]][1] = v[buckets[i][k]][1] + newTimeStepSize * force1[buckets[i][k]] / mass[buckets[i][k]];
                        v[buckets[i][k]][2] = v[buckets[i][k]][2] + newTimeStepSize * force2[buckets[i][k]] / mass[buckets[i][k]];
                }

            }

        }
    }

    // If there were collisions merge the collided bodies now
    if (collisionTotal > 0) {
        for (int ic = 0; ic < collisionTotal; ic++) {
            int first = collision1[ic];
            int second = collision2[ic];

            x[first][0] = (x[first][0] + x[second][0]) / 2;
            x[first][1] = (x[first][1] + x[second][1]) / 2;
            x[first][2] = (x[first][2] + x[second][2]) / 2;

            v[first][0] = ((mass[first] / (mass[first] + mass[second])) * v[first][0]) + ((mass[second] / (mass[first] + mass[second])) * v[second][0]);
            v[first][1] = ((mass[first] / (mass[first] + mass[second])) * v[first][1]) + ((mass[second] / (mass[first] + mass[second])) * v[second][1]);
            v[first][2] = ((mass[first] / (mass[first] + mass[second])) * v[first][2]) + ((mass[second] / (mass[first] + mass[second])) * v[second][2]);
            mass[first] = mass[first] + mass[second];

            NumberOfBodies -= 1;
            // Shift array to delete the second colliding object as fused object has already been stored in the first
            // This cannot be done in the main loop while other threads are using the prorerties of the bodies
            for (int m = second; m < NumberOfBodies; m++) {
                x[m] = x[m + 1]; 
                v[m] = v[m + 1];
                mass[m] = mass[m + 1];
            }
            if (NumberOfBodies == 1) {
                tFinal = t; //print position of last object
                std::cout << x[0][0]
                << ", " << x[0][1]
                << ", " << x[0][2]
                << std::endl;
            }
        }
    }

    t = t + timeStepSize;

    delete[] force0;
    delete[] force1;
    delete[] force2;
    delete[] forces;
    delete[] collision1;
    delete[] collision2;
    delete[] velocities;
    for (int i = 0; i < numBuckets; i++)
        delete[] buckets[i];
    delete[] buckets;
}

/**
 * Main routine.
 *
 * Not to be changed in assignment.
 */
int main(int argc, char** argv) {
    if (argc == 1) {
        std::cerr << "usage: " + std::string(argv[0]) + " snapshot final-time dt objects" << std::endl
            << "  snapshot        interval after how many time units to plot. Use 0 to switch off plotting" << std::endl
            << "  final-time      simulated time (greater 0)" << std::endl
            << "  dt              time step size (greater 0)" << std::endl
            << std::endl
            << "Examples:" << std::endl
            << "0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0 \t One body moving form the coordinate system's centre along x axis with speed 1" << std::endl
            << "0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0     0.0 1.0 0.0  1.0 0.0 0.0  1.0  \t One spiralling around the other one" << std::endl
            << "0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0 \t Three body setup from first lecture" << std::endl
            << "0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0     2.0 1.0 0.0  0.0 0.0 0.0  1.0     2.0 0.0 1.0  0.0 0.0 0.0  1.0 \t Five body setup" << std::endl
            << std::endl
            << "In this naive code, only the first body moves" << std::endl;

        return -1;
    }
    else if ((argc - 4) % 7 != 0) {
        std::cerr << "error in arguments: each planet is given by seven entries (position, velocity, mass)" << std::endl;
        std::cerr << "got " << argc << " arguments (three of them are reserved)" << std::endl;
        std::cerr << "run without arguments for usage instruction" << std::endl;
        return -2;
    }

    std::cout << std::setprecision(15);

    setUp(argc, argv);

    openParaviewVideoFile();

    int snapshotCounter = 0;
    if (t > tPlot) {
        printParaviewSnapshot();
        std::cout << "plotted initial setup" << std::endl;
        tPlot = tPlotDelta;
    }

    int timeStepCounter = 0;
    while (t <= tFinal) {
        updateBody();
        timeStepCounter++;
        if (t >= tPlot) {
            printParaviewSnapshot();
            std::cout << "plot next snapshot"
                << ",\t time step=" << timeStepCounter
                << ",\t t=" << t
                << ",\t dt=" << timeStepSize
                << ",\t v_max=" << maxV
                << ",\t dx_min=" << minDx
                << std::endl;

            tPlot += tPlotDelta;
        }
    }

    closeParaviewVideoFile();

    std::cout << "exit normally" << std::endl;

    return 0;
}