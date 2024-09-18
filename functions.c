#include <stdio.h>          // standard c input/output library
#include <stdlib.h>         // standard c library
#include <string.h>
#include <math.h>
#include "functions.h"

void initBody(struct Body *b, char *name, double mass, double *r, double *v) {
    strncpy(b->name, name, 12); // Copy name to struct
    b->m = mass;                // Set mass of the body
    
    for (int i = 0; i < 3; i++) {
        b->r[i] = r[i];         // Assign each element of the position vector
        b->v[i] = v[i];         // Assign each element of the velocity vector
        b->a[i] = 0;            // Initialize acceleration to 0
        b->pa[i] = 0;           // Initialize previous acceleration to 0
    }
}

void calculateAccelerations(struct Body *b) {
    double x_ij, y_ij, z_ij, r_ij2, r_ij, a_ij, a_ji;
    Utot = 0;   // Reset potential energy
    // Reset accelerations
    for (int i = 0; i < numberOfBodies; i++) {
        b[i].a[0] = 0;
        b[i].a[1] = 0;
        b[i].a[2] = 0;
    }
    // Loop over the pair interactions
    for (int i = 0; i < numberOfBodies; i++) {
        for (int j = i + 1; j < numberOfBodies; j++) {
            // Calculate the distance between the bodies
            x_ij = b[i].r[0] - b[j].r[0];
            y_ij = b[i].r[1] - b[j].r[1];
            z_ij = b[i].r[2] - b[j].r[2];
            r_ij2 = x_ij * x_ij + y_ij * y_ij + z_ij * z_ij;
            r_ij = sqrt(r_ij2);

            // Calculate the common factor of acceleration
            a_ij = -G / r_ij2;

            // Calculate the acceleration components for body i and j
            b[i].a[0] += a_ij * (x_ij / r_ij) * b[j].m;
            b[i].a[1] += a_ij * (y_ij / r_ij) * b[j].m;
            b[i].a[2] += a_ij * (z_ij / r_ij) * b[j].m;

            b[j].a[0] -= a_ij * (x_ij / r_ij) * b[i].m;
            b[j].a[1] -= a_ij * (y_ij / r_ij) * b[i].m;
            b[j].a[2] -= a_ij * (z_ij / r_ij) * b[i].m;

            // Calculate the potential energy
            Utot += -G * b[i].m * b[j].m / r_ij;
        }
    }
}

void velocityVerlet(struct Body *b) {
    // Update positions and store previous accelerations
    for(int i = 0; i < numberOfBodies; i++) {
        b[i].r[0] += b[i].v[0] * dt + half_dt2 * b[i].a[0];
        b[i].r[1] += b[i].v[1] * dt + half_dt2 * b[i].a[1];
        b[i].r[2] += b[i].v[2] * dt + half_dt2 * b[i].a[2];

        b[i].pa[0] = b[i].a[0];
        b[i].pa[1] = b[i].a[1];
        b[i].pa[2] = b[i].a[2];
    }

    // Calculate new accelerations
    calculateAccelerations(b);

    // Update velocities
    for(int i = 0; i < numberOfBodies; i++) {
        b[i].v[0] += half_dt * (b[i].pa[0] + b[i].a[0]);
        b[i].v[1] += half_dt * (b[i].pa[1] + b[i].a[1]);
        b[i].v[2] += half_dt * (b[i].pa[2] + b[i].a[2]);
    }    
}

void writePositions(struct Body *b, FILE *posFile) {
    for (int i = 0; i < numberOfBodies; i++) {
        fprintf(posFile, "%-5d %-12s %16.5g %16.5g %16.5g\n", t, b[i].name, b[i].r[0], b[i].r[1], b[i].r[2]);
    }
}

void calculateEnergy(struct Body *b, FILE *energyFile) {
    double v2;
    double Ktot = 0;    // Reset total kinetic energy

    for(int i = 0; i < numberOfBodies; i++) {
        v2 = b[i].v[0] * b[i].v[0] + b[i].v[1] * b[i].v[1] + b[i].v[2] * b[i].v[2]; // Calculate velocity squared
        Ktot += 0.5 * b[i].m * v2; // Calculate kinetic energy
    }

    double Etot = Ktot + Utot; // Calculate total energy
    writeEnergies(Ktot, Utot, Etot, t, energyFile); // Write energies to file
}

void initEnergyFile(FILE *energyFile) {
    fprintf(energyFile, "%-5s\t%-15s\t%-15s\t%-15s\n", "t", "Ktot", "Utot", "Etot");
}

void writeEnergies(double Ktot, double Utot, double Etot, int t, FILE *energyFile) {
    fprintf(energyFile, "%-5d\t%-15.10e\t%-15.10e\t%-15.10e\n", t, Ktot, Utot, Etot);
}

void readData(const char *filename, struct Body *b) {
    FILE *file = fopen(filename, "r");     // Pointer to file location
    const int n = 1024;                    // Maximum number of characters per entry
    char line[n];                          // Pointer to array when read string is stored

    // Specify data types of each column
    int id;
    char name[20];
    double mass,rx,ry,rz,vx,vy,vz;

    // Error if the file cannot be found
    if(file == NULL) {
        printf("Error opening file.\n");
        exit(1);
    }

    int bodyCounter = 0;
    // // Reads every row of data from file
    while (fgets(line, n, file) != NULL) {
        // Ignore first line which starts with hashtag
        if (line[0] != '#') {
            // Read raw string of data and stores in provided variables
            sscanf(line, "%d %s %lf %lf %lf %lf %lf %lf %lf", &id, name, &mass, &rx, &ry, &rz, &vx, &vy, &vz); 

            // Make array with position and vector components
            double r[3] = {rx, ry, rz};
            double v[3] = {vx, vy, vz};

            // Use initialize body function to assign data to body struct
            initBody(&b[bodyCounter], name, mass, r, v);
            
            // Next body
            bodyCounter++;
        }
    }

    fclose(file);
}

void calculateDeviation(struct Body *simulated, struct Body *real, FILE *deviationFile) {
    double deviation[3];
    fprintf(deviationFile, "Name,X,Y,Z\n");
    for (int i = 0; i < numberOfBodies; i++) {
        deviation[0] = fabs(simulated[i].r[0] - real[i].r[0]) / fabs(real[i].r[0]) * 100;
        deviation[1] = fabs(simulated[i].r[1] - real[i].r[1]) / fabs(real[i].r[1]) * 100;
        deviation[2] = fabs(simulated[i].r[2] - real[i].r[2]) / fabs(real[i].r[2]) * 100;

        fprintf(deviationFile, "%s,%.5e,%.5e,%.5e\n", simulated[i].name, deviation[0], deviation[1], deviation[2]);
    }
}

void symplecticEuler(struct Body *b) {
    // Update velocity
    for (int i = 0; i < numberOfBodies; i++) {
        b[i].v[0] += b[i].a[0] * dt;
        b[i].v[1] += b[i].a[1] * dt;
        b[i].v[2] += b[i].a[2] * dt;
    }

    // Update position
    for (int i = 0; i < numberOfBodies; i++) {
        b[i].r[0] += b[i].v[0] * dt;
        b[i].r[1] += b[i].v[1] * dt;
        b[i].r[2] += b[i].v[2] * dt;
    }

    calculateAccelerations(b); // Update acceleration
}

void explicitEuler(struct Body *b) {
    // Update positions
    for (int i = 0; i < numberOfBodies; i++) {
        b[i].r[0] += b[i].v[0] * dt;
        b[i].r[1] += b[i].v[1] * dt;
        b[i].r[2] += b[i].v[2] * dt;
    }

    // Update velocities
    for (int i = 0; i < numberOfBodies; i++) {
        b[i].v[0] += b[i].a[0] * dt;
        b[i].v[1] += b[i].a[1] * dt;
        b[i].v[2] += b[i].a[2] * dt;
    }

    calculateAccelerations(b); // Update acceleration
}