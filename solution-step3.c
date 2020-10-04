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


double t          = 0;
double tFinal     = 0;
double tPlot      = 0;
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
double*  mass;

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
  NumberOfBodies = (argc-4) / 7;

  x    = new double*[NumberOfBodies];
  v    = new double*[NumberOfBodies];
  mass = new double [NumberOfBodies];

  int readArgument = 1;

  tPlotDelta   = std::stof(argv[readArgument]); readArgument++;
  tFinal       = std::stof(argv[readArgument]); readArgument++;
  timeStepSize = std::stof(argv[readArgument]); readArgument++;

  for (int i=0; i<NumberOfBodies; i++) {
    x[i] = new double[3];
    v[i] = new double[3];

    x[i][0] = std::stof(argv[readArgument]); readArgument++;
    x[i][1] = std::stof(argv[readArgument]); readArgument++;
    x[i][2] = std::stof(argv[readArgument]); readArgument++;

    v[i][0] = std::stof(argv[readArgument]); readArgument++;
    v[i][1] = std::stof(argv[readArgument]); readArgument++;
    v[i][2] = std::stof(argv[readArgument]); readArgument++;

    mass[i] = std::stof(argv[readArgument]); readArgument++;

    if (mass[i]<=0.0 ) {
      std::cerr << "invalid mass for body " << i << std::endl;
      exit(-2);
    }
  }

  std::cout << "created setup with " << NumberOfBodies << " bodies" << std::endl;
  
  if (tPlotDelta<=0.0) {
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
  videoFile.open( "result.pvd" );
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
  filename << "result-" << counter <<  ".vtp";
  std::ofstream out( filename.str().c_str() );
  out << "<VTKFile type=\"PolyData\" >" << std::endl
      << "<PolyData>" << std::endl
      << " <Piece NumberOfPoints=\"" << NumberOfBodies << "\">" << std::endl
      << "  <Points>" << std::endl
      << "   <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">";
//      << "   <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">";

  for (int i=0; i<NumberOfBodies; i++) {
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
      << "</VTKFile>"  << std::endl;

  videoFile << "<DataSet timestep=\"" << counter << "\" group=\"\" part=\"0\" file=\"" << filename.str() << "\"/>" << std::endl;
}



/**
 * This is the only operation you are allowed to change in the assignment.
 */
void updateBody() {
  maxV   = 0.0;
  minDx  = std::numeric_limits<double>::max();

  // force0 = force along x direction
  // force1 = force along y direction
  // force2 = force along z direction
  double* force0 = new double[NumberOfBodies];
  double* force1 = new double[NumberOfBodies];
  double* force2 = new double[NumberOfBodies];
  double* forces = new double[NumberOfBodies*3];
  int* collisions = new int[NumberOfBodies];
  double* velocities = new double[NumberOfBodies];
  int numBuckets = 10;
  int * * buckets = new int *[numBuckets];
  for (int i=0;i<numBuckets;++i){
    buckets[i] = new int[NumberOfBodies+1];
  }
  double vbucket = 0.0;
  bool assigned = false;
  float newTimeStepSize = timeStepSize;

// first element in each bucket counts the number of objects in the bucket
 for (int i=0;i<numBuckets;i++){
   buckets[i][0] = 0.0;
 }

  // storing the velocities and max velocity of the particles
  for (int i=0;i<NumberOfBodies;i++){
    velocities[i] = sqrt( v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2] );
    if (i==0){
      maxV = std::sqrt( v[0][0]*v[0][0] + v[0][1]*v[0][1] + v[0][2]*v[0][2] );
    } else {
      maxV = std::max(maxV,sqrt( v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2] ));
    }
  }
  
  vbucket = maxV/numBuckets;
  for (int i=0; i<NumberOfBodies;i++){
      assigned = false;
        for (int j=0; j<numBuckets;j++){
          //adding an element to the bucket if its velocity is bounded by the bucket
          if ((j)*vbucket <= velocities[i] && velocities[i]<((j+1)*vbucket)){
            buckets[j][int(buckets[j][0])+1] = i;
            buckets[j][0]+=1;//increment the counter of objects in the bucket
            assigned = true;
            break;
          }
        }
      if (assigned == false) { //must be maxV and in last bucket
        buckets[numBuckets-1][int(buckets[numBuckets-1][0])+1] = i;
        buckets[numBuckets-1][0]+=1;
      }
  }
  //loop through the buckets
  for (int i=0; i<numBuckets; i++){
    if (i>0){ //halve the time step size each iteration
      newTimeStepSize = newTimeStepSize/2;
    }
 
    for (int j=0;j<pow(2,i);j++){ //number of times to run for each bucket is 2**i

      for (int k=1; k<buckets[i][0]+1; k++) { //for all objects in bucket i

        collisions[buckets[i][k]] = -1;
        force0[buckets[i][k]] = 0.0;
        force1[buckets[i][k]] = 0.0;
        force2[buckets[i][k]] = 0.0;
        

        for (int l=0; l<NumberOfBodies; l++) {

            if (buckets[i][k]!=l) {
                const double distance = sqrt(
                (x[buckets[i][k]][0]-x[l][0]) * (x[buckets[i][k]][0]-x[l][0]) +
                (x[buckets[i][k]][1]-x[l][1]) * (x[buckets[i][k]][1]-x[l][1]) +
                (x[buckets[i][k]][2]-x[l][2]) * (x[buckets[i][k]][2]-x[l][2])
                );

                // x,y,z forces acting on particle 0
                force0[buckets[i][k]] += (x[l][0]-x[buckets[i][k]][0]) * mass[l]*mass[buckets[i][k]] / distance / distance / distance ;
                force1[buckets[i][k]] += (x[l][1]-x[buckets[i][k]][1]) * mass[l]*mass[buckets[i][k]] / distance / distance / distance ;
                force2[buckets[i][k]] += (x[l][2]-x[buckets[i][k]][2]) * mass[l]*mass[buckets[i][k]] / distance / distance / distance ;

                if (k<l) {// no need to calculate twice
                    minDx = std::min( minDx,distance );
                    //collision if distance between 2 objects less than sum of radii
                    if (distance <= 0.01) {
                      collisions[buckets[i][k]] = l; //store position of object collided with
                    }
                }            
            }
        }
        //storing the forces acting on each object
        forces[3*buckets[i][k]] = force0[buckets[i][k]];
        forces[(3*buckets[i][k])+1] = force1[buckets[i][k]];
        forces[(3*buckets[i][k])+2] = force2[buckets[i][k]];
      }
      //updating the positions, velocities, maxV and mass depending on if the particle collided
      for (int k=1; k<buckets[i][0]+1; k++){
        if (collisions[buckets[i][k]] == -1){ //no collision - update as normal
            x[buckets[i][k]][0] = x[buckets[i][k]][0] + newTimeStepSize * v[buckets[i][k]][0];
            x[buckets[i][k]][1] = x[buckets[i][k]][1] + newTimeStepSize * v[buckets[i][k]][1];
            x[buckets[i][k]][2] = x[buckets[i][k]][2] + newTimeStepSize * v[buckets[i][k]][2];

            v[buckets[i][k]][0] = v[buckets[i][k]][0] + newTimeStepSize * force0[buckets[i][k]] / mass[buckets[i][k]];
            v[buckets[i][k]][1] = v[buckets[i][k]][1] + newTimeStepSize * force1[buckets[i][k]] / mass[buckets[i][k]];
            v[buckets[i][k]][2] = v[buckets[i][k]][2] + newTimeStepSize * force2[buckets[i][k]] / mass[buckets[i][k]];

        } else { //collision - calculate parameters of fused object
          x[buckets[i][k]][0] = (x[buckets[i][k]][0] + x[collisions[buckets[i][k]]][0])/2;
          x[buckets[i][k]][1] = (x[buckets[i][k]][1] + x[collisions[buckets[i][k]]][1])/2;
          x[buckets[i][k]][2] = (x[buckets[i][k]][2] + x[collisions[buckets[i][k]]][2])/2;
          v[buckets[i][k]][0] = ((mass[buckets[i][k]]/(mass[buckets[i][k]]+mass[collisions[buckets[i][k]]]))*v[buckets[i][k]][0]) + ((mass[collisions[buckets[i][k]]]/(mass[buckets[i][k]]+mass[collisions[buckets[i][k]]]))*v[collisions[buckets[i][k]]][0]);
          v[buckets[i][k]][1] = ((mass[buckets[i][k]]/(mass[buckets[i][k]]+mass[collisions[buckets[i][k]]]))*v[buckets[i][k]][1]) + ((mass[collisions[buckets[i][k]]]/(mass[buckets[i][k]]+mass[collisions[buckets[i][k]]]))*v[collisions[buckets[i][k]]][1]);
          v[buckets[i][k]][2] = ((mass[buckets[i][k]]/(mass[buckets[i][k]]+mass[collisions[buckets[i][k]]]))*v[buckets[i][k]][2]) + ((mass[collisions[buckets[i][k]]]/(mass[buckets[i][k]]+mass[collisions[buckets[i][k]]]))*v[collisions[buckets[i][k]]][2]);
          mass[buckets[i][k]] = mass[buckets[i][k]] + mass[collisions[buckets[i][k]]];
          NumberOfBodies -= 1;
          for (int m = collisions[buckets[i][k]];m<NumberOfBodies;m++) {
            x[m] = x[m+1]; //shift array to delete one of the colliding objects
            v[m] = v[m+1]; //as fused object has already been stored in the other
            mass[m] = mass[m+1];
          }
          if (NumberOfBodies == 1) {
            tFinal = t; //print out position of last object and terminate program
            std::cout << x[0][0]
                  << ", " << x[0][1]
                  << ", " << x[0][2]
                  << std::endl;
            }
            // std::cout << x[0][0] //to print out all collisions
            //     << ", " << x[0][1]
            //     << ", " << x[0][2]
            //      << std::endl;
        }
        maxV = std::max(maxV,sqrt( v[buckets[i][k]][0]*v[buckets[i][k]][0] + v[buckets[i][k]][1]*v[buckets[i][k]][1] + v[buckets[i][k]][2]*v[buckets[i][k]][2] ));
      }
    }
  }

  t += timeStepSize;

  delete[] force0;
  delete[] force1;
  delete[] force2;
  delete[] forces;
  delete[] collisions;
  delete[] velocities;
  for(int i = 0; i < numBuckets; i++)
    delete[] buckets[i];
  delete[] buckets;
}


/**
 * Main routine.
 *
 * Not to be changed in assignment.
 */
int main(int argc, char** argv) {
  if (argc==1) {
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
  else if ( (argc-4)%7!=0 ) {
    std::cerr << "error in arguments: each planet is given by seven entries (position, velocity, mass)" << std::endl;
    std::cerr << "got " << argc << " arguments (three of them are reserved)" << std::endl;
    std::cerr << "run without arguments for usage instruction" << std::endl;
    return -2;
  }

  std::cout << std::setprecision(15);

  setUp(argc,argv);

  openParaviewVideoFile();

  int snapshotCounter = 0;
  if (t > tPlot) {
    printParaviewSnapshot();
    std::cout << "plotted initial setup" << std::endl;
    tPlot = tPlotDelta;
  }

  int timeStepCounter = 0;
  while (t<=tFinal) {
    updateBody();
    timeStepCounter++;
    if (t >= tPlot) {
      printParaviewSnapshot();
      std::cout << "plot next snapshot"
    		    << ",\t time step=" << timeStepCounter
    		    << ",\t t="         << t
				<< ",\t dt="        << timeStepSize
				<< ",\t v_max="     << maxV
				<< ",\t dx_min="    << minDx
				<< std::endl;

      tPlot += tPlotDelta;
    }
  }

  closeParaviewVideoFile();

  return 0;
}