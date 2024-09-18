#ifndef FUNCTIONS_H
#define FUNCTIONS_H

// Long-term simulation parameters
// #define Tmax (200000LL * 24 * 3600LL)
// #define numberOfBodies 10
// #define dt (100LL * 24 * 3600LL)

// Simulation parameters
#define Tmax 365 * 24 * 3600    // One year in seconds
#define numberOfBodies 230      // Number of bodies in the simulation
#define dt 3600                 // Time step in seconds
#define G 6.67430E-11           // Gravitational constant

extern int steps;
extern int t;
extern const double half_dt;
extern const double half_dt2;
extern double Utot;

// Define the Body structure
struct Body {
    char name[20];
    double m;       // Mass
    double r[3];    // Position
    double v[3];    // Velocity
    double a[3];    // Acceleration
    double pa[3];   // Previous acceleration
};

// Function declarations
void readData(const char *filename, struct Body *b);
void initBody(struct Body *b, char *name, double mass, double *r, double *v);

void velocityVerlet(struct Body *b);
void calculateAccelerations(struct Body *b);

void writePositions(struct Body *b, FILE *posFile);
void calculateEnergy(struct Body *b, FILE *energyFile);
void initEnergyFile(FILE *energyFile);
void writeEnergies(double Ktot, double Utot, double Etot, int t, FILE *energyFile);
void calculateDeviation(struct Body *simulated, struct Body *real, FILE *deviationFile);

void symplecticEuler(struct Body *b);
void explicitEuler(struct Body *b);

#endif // FUNCTIONS_H