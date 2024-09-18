#include <stdio.h>          // standard c input/output library
#include <stdlib.h>         // standard c library
#include <string.h>
#include <math.h>

#include "functions.h"

int steps = Tmax / dt;
int t = 0;
const double half_dt = 0.5 * dt;
const double half_dt2 = 0.5 * dt * dt;
double Utot = 0;

int main(void) {
    struct Body b[numberOfBodies];

    readData("solar_system_13sept2021.dat", b);

    FILE *posFile = fopen("positions.dat", "w");
    FILE *energyFile = fopen("energies.dat", "w");

    initEnergyFile(energyFile);
    calculateAccelerations(b);

    // Ask user to choose an algorithm
    int choice;
    printf("Choose an integration algorithm:\n");
    printf("1. Velocity Verlet\n");
    printf("2. Symplectic Euler\n");
    printf("3. Explicit Euler\n");
    printf("Enter your choice (1-3): ");
    scanf("%d", &choice);

    for (t = 0; t < steps; t++) {
        // Call the selected integration function
        switch (choice) {
            case 1:
                velocityVerlet(b);
                break;
            case 2:
                symplecticEuler(b);
                break;
            case 3:
                explicitEuler(b);
                break;
            default:
                fprintf(stderr, "Invalid choice. Exiting.\n");
                fclose(posFile);
                fclose(energyFile);
                return 1;
        }

        writePositions(b, posFile);  // Write positions to file
        calculateEnergy(b, energyFile); // Calculate and write energies
    }

    fclose(posFile);
    fclose(energyFile);

    struct Body real[numberOfBodies];
    FILE *deviationFile = fopen("deviation.dat", "w");
    readData("solar_system_13sept2022.dat", real);
    calculateDeviation(b, real, deviationFile);
    
	return (0);
}