/********************************************************************************************************************

Ziad RAMMAL
Xing WEI
Jinhao XU

Projet niveau A : Genetique

********************************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

/********************************************************************************************************************/

// Initialisation des constantes :

#define Nb_bits 8      // Nombre de bits
#define Nb_lambda 1024 // Nombre de lambda

const int POP_SIZE = 1000;      // Nombre de la population
const double MUTATION_RATE = 10; // Pourcentage de mutation
const double pow2_Nb_bits_minus_1 = (1 << Nb_bits) - 1;
//const double pow2_Nb_bits_minus_1 = (pow(2, Nb_bits) - 1);

// Initialisation des bornes des variables :

const double Y0_min = 0.21;
const double Y0_max = 0.23;
const double B_min = 12.5;
const double B_max = 13.5;
const double lambda0_min = 6561.0;
const double lambda0_max = 6563.0;
const double deltaLambdaL_min = 9.7;
const double deltaLambdaL_max = 10.2;

/********************************************************************************************************************/

// Définir des type de variables

typedef double tab[Nb_lambda];

typedef struct
{
    long int code_bin_Y0;
    long int code_bin_B;
    long int code_bin_lambda0;
    long int code_bin_deltaLambdaL;
    double err;
} Individue_bin;

typedef struct
{
    double Y0;
    double B;
    double lambda0;
    double deltaLambdaL;
    double err;
} Individue_dec;


/********************************************************************************************************************/

// Appel des fonctions :

void readData(tab lambda, tab bruit, tab lisse, int *Np_lamda); // Fonction qui lit les données du fichier "Data_Genetique.txt"

long convert_to_decimal(long num);

long convert_to_binaire(long int num);

long convert_to_binaire_optimized(long int num);

double Y0_en_dec(int Y0_bin);

double B_en_dec(int B_bin);

double lambda0_en_dec(int lambda0_bin);

double deltaLambdaL_en_dec(int deltaLambdaL_bin);

void F_lorentzienne(Individue_bin ind_bin, const tab lambda, tab J);

void tracee_finale(double *lambda, double *J);

void calc_err(Individue_bin *ind_bin, const tab bruit, const tab J);

void Ind_bin_to_Ind_dec(Individue_bin *ind_bin, Individue_dec *ind_dec);

void generer_individus(Individue_dec *individus, Individue_bin *individus_bin);

long int Y0_en_bin(double Y0_dec);

long int B_en_bin(double B_dec);

long int lambda0_en_bin(double lambda0_dec);

long int deltaLambdaL_en_bin(double deltaLambdaL_dec);

void Ind_dec_to_Ind_bin(Individue_bin *ind_bin, Individue_dec *ind_dec);

int *generer_nombres(int taille);

void Duel(Individue_bin individues[], Individue_bin nouvellePop[], int taille, int *new_pop_size);

int exposant(int base, int puissance);

void mutate_chromosome(long int *individu);

void mutate(Individue_bin *individu);

void croisement_chromosome(long int ind1, long int ind2, long int *enfant1, long int *enfant2);

void croisement(Individue_bin parent1, Individue_bin parent2, Individue_bin *enfant1, Individue_bin *enfant2);

double densite(double deltaLambdaL);

/********************************************************************************************************************/

// Fonction principale :

int main(void)
{
    int N, new_pop_size = 0;
    double densite_f;
    tab lambda, bruit, lisse, J;
    Individue_bin ind_bin[POP_SIZE];
    Individue_dec individus[POP_SIZE];
    Individue_dec best_dec;
    Individue_bin enfant1, enfant2;

    // Set seed for random number generation
    srand(time(NULL));

    // Lecture des données du fichier
    readData(lambda, bruit, lisse, &N);

    // Génération des individus
    generer_individus(individus, ind_bin);

    int temp_POP_SIZE = POP_SIZE;

    // Evolution de la population
    //while (temp_POP_SIZE > 1)
    for (int i = 0; i < 100; i++)
    {
        for (int j = 0; j < temp_POP_SIZE; j++)
     {
        // Calculate Lorentzian function for each individual
         F_lorentzienne(ind_bin[j], lambda, J);
         // Calculate error for each individual
         calc_err(&ind_bin[j] , bruit, J);
     }
        // Sélection
        Individue_bin nouvellePop[temp_POP_SIZE/2];
        //int *nombres = generer_nombres(temp_POP_SIZE);
        Duel(ind_bin, nouvellePop, temp_POP_SIZE, &new_pop_size);

        // Croisement
        for (int j = 0; j <new_pop_size; j+=2)
        {
            croisement(nouvellePop[j], nouvellePop[j + 1], &enfant1, &enfant2);

            // Mutation
            mutate(&enfant1);
            mutate(&enfant2);

            // Ajout des enfants mutés à la nouvelle population
            nouvellePop[j] = enfant1;
            nouvellePop[j + 1] = enfant2;
        }

        // Copie de la nouvelle population dans la population actuelle
        memcpy(ind_bin, nouvellePop, new_pop_size * sizeof(Individue_bin));
        temp_POP_SIZE = new_pop_size;
    }

    // Calcul de l'erreur de l'individu le plus performant
    Individue_bin best_ind = ind_bin[0];
    for (int i = 1; i < new_pop_size; i++)
    {
        if (ind_bin[i].err < best_ind.err)
        {
            best_ind = ind_bin[i];
        }
    }
    Ind_bin_to_Ind_dec(&best_ind,&best_dec);
    // Calcul de la fonction de Lorenzienne et traçage de la courbe
    F_lorentzienne(best_ind, lambda, J);
    tracee_finale(lambda, J);

    densite_f = densite (best_dec.deltaLambdaL);

    printf("L'algorithme génétique binaire est en cours d'exécution...\n");
    printf("Individu : Y0= %f\t ; B = %f\t ; lambda0 = %f\t ; deltaL = %f\n", best_dec.Y0, best_dec.B, best_dec.lambda0, best_dec.deltaLambdaL);
    printf("densite = %lf\n",densite_f);

    return 0;
}

/********************************************************************************************************************/

// Initialisation des fonctions :

void readData(tab lambda, tab bruit, tab lisse, int *Np_lamda) // Fonction qui lit les données du fichier "Data_Genetique.txt"
{
    int i = 0;
    FILE *file = fopen("Data_Genetique.txt", "r"); // Ouverture du fichier
    if (file == NULL)                              // Vérifier que le fichier n'est pas vide et bien ouvert
    {
        printf("Erreur lors de l'ouverture du fichier\n");
        exit(1);
    }

    while (!feof(file)) // Lecture des données jusqu'à la fin du dossier
    {
        fscanf(file, "%lf %lf %lf\n", &lambda[i], &bruit[i], &lisse[i]);
        i++; // Stocker les valeures dans des tableaux
    }
    fclose(file); // Fermeture du fichier
    *Np_lamda = i;
}
/*________________________________________*/

long convert_to_decimal(long num)
{
    int dec = 0, rem = 0, place = 0;
    while (num)
    {
        rem = num % 10;
        if (rem) { // Check if remainder exists
            dec += exposant(2, place);
        }

        num = num / 10;
        place++;
    }
    return dec;
}
/*________________________________________*/

long convert_to_binaire(long int num)
{
    int rem = 0, place = 1; // Remainder, Place value variables
    int bin = 0; // Binary result variable

    while (num) // If num is not equal to 0
    {
        rem = num % 2; // Find remainder when num is divided by 2
        num = num / 2; // Divide num by 2 and update
        bin = bin + (rem * place); // Update binary result
        place = place * 10; // Update place value for next iteration
    }
    return bin; // Return binary result
}
/*________________________________________*/

double Y0_en_dec(int Y0_bin)
{
    // Convert binary to decimal
    int Y0_bin_dec = convert_to_decimal(Y0_bin);

    // Decode binary value to physical quantity
    //double Y0_dec = (((Y0_max - Y0_min) / pow2_Nb_bits_minus_1) * Y0_bin_dec) + Y0_min;
    double Y0_dec = ((((Y0_max - Y0_min) / pow2_Nb_bits_minus_1) * Y0_bin_dec) + Y0_min);

    return Y0_dec;
}
/*________________________________________*/

double B_en_dec(int B_bin)
{
    // Convert binary to decimal
    int B_bin_dec = convert_to_decimal(B_bin);

    // Perform backward binary to decimal mapping
    return ((((B_max - B_min) / (pow2_Nb_bits_minus_1) * B_bin_dec) + B_min));
}
/*________________________________________*/

double lambda0_en_dec(int lambda0_bin)
{
    // Convert binary to decimal
    int lambda0_bin_dec = convert_to_decimal(lambda0_bin);

    // Encode lambda0 using formula
    double lambda0_en = ((((lambda0_max - lambda0_min) / (pow2_Nb_bits_minus_1) * lambda0_bin_dec) + lambda0_min));

    return lambda0_en;
}
/*________________________________________*/

double deltaLambdaL_en_dec(int deltaLambdaL_bin)
{
    // Convert binary to decimal
    int deltaLambdaL_bin_dec = convert_to_decimal(deltaLambdaL_bin);

    // Convert bin to deltaLambdaL value
    return ((((deltaLambdaL_max - deltaLambdaL_min) / (pow2_Nb_bits_minus_1) * deltaLambdaL_bin_dec) + deltaLambdaL_min));
}
/*________________________________________*/

void F_lorentzienne(Individue_bin ind_bin, const tab lambda, tab J)
{
    // Convert binary to decimal values
    double Y0 = Y0_en_dec(ind_bin.code_bin_Y0);
    double B = B_en_dec(ind_bin.code_bin_B);
    double lambda0 = lambda0_en_dec(ind_bin.code_bin_lambda0);
    double deltaLambdaL = deltaLambdaL_en_dec(ind_bin.code_bin_deltaLambdaL);

    // Compute intensity at each wavelength
    for (int i = 0; i < Nb_lambda; i++)
    {
        J[i] = Y0 + (B / (((lambda[i] - lambda0) * (lambda[i] - lambda0)) + (((deltaLambdaL) * (deltaLambdaL)) / 4)));
    }
}
/*________________________________________*/

void tracee_finale(double *lambda, double *J)
{
    FILE *file;
    file = fopen("Tracee_Finale.txt", "w"); // Ouverture du fichier en écriture
    if (file == NULL)
    {
        printf("Erreur lors de la création du fichier\n"); // Gestion de l'erreur
        return;
    }
    for (int i = 0; i < Nb_lambda; i++)
    {
        fprintf(file, "%lf %lf\n", lambda[i], J[i]); // Ecriture des valeurs dans le fichier
    }
    fclose(file); // Fermeture du fichier
}
/*________________________________________*/

void calc_err(Individue_bin *ind_bin, const tab bruit, const tab J)
{
    // Initialize error sum to 0
    double sum = 0.0;

    // Loop through all wavelengths
    for (int i = 0; i < Nb_lambda; i++)
    {
        // Calculate squared difference and add to sum
        sum = sum + ((bruit[i] - J[i]) * (bruit[i] - J[i]));
    }

    // Store total error in individual's error attribute
    ind_bin->err = sum;
}
/*________________________________________*/

void Ind_bin_to_Ind_dec(Individue_bin *ind_bin, Individue_dec *ind_dec)
{
    // Convert binary to decimal Y0
    ind_dec->Y0 = Y0_en_dec(ind_bin->code_bin_Y0);
    
    // Convert binary to decimal B
    ind_dec->B = B_en_dec(ind_bin->code_bin_B);
    
    // Convert binary to decimal lambda0
    ind_dec->lambda0 = lambda0_en_dec(ind_bin->code_bin_lambda0);
    
    // Convert binary to decimal deltaLambdaL
    ind_dec->deltaLambdaL = deltaLambdaL_en_dec(ind_bin->code_bin_deltaLambdaL);
    
    // Copy error from binary to decimal
    ind_dec->err = ind_bin->err;
}
/*________________________________________*/

void generer_individus(Individue_dec *individus, Individue_bin *individus_bin)
{
    // Set seed for random number generation

    // Generate population size number of individuals
    for (int i = 0; i < POP_SIZE; i++)
    {
        // Generate random value Y0 within range
        int Y0 = rand() % (int)(pow2_Nb_bits_minus_1);
        individus_bin[i].code_bin_Y0 = convert_to_binaire(Y0);
        individus[i].Y0 = ((((Y0_max - Y0_min) / (pow2_Nb_bits_minus_1) * Y0) + Y0_min));

        // Generate random value B within range
        int B = rand() % (int)(pow2_Nb_bits_minus_1);
        individus_bin[i].code_bin_B = convert_to_binaire(B);
        individus[i].B = ((((B_max - B_min) / (pow2_Nb_bits_minus_1) * B) + B_min));

        // Generate random value lambda0 within range
        int lambda0 = rand() % (int)(pow2_Nb_bits_minus_1);
        individus_bin[i].code_bin_lambda0 = convert_to_binaire(lambda0);
        individus[i].lambda0 = ((((lambda0_max - lambda0_min) / (pow2_Nb_bits_minus_1) * lambda0) + lambda0_min));

        // Generate random value deltaLambdaL within range
        int deltaLambdaL = rand() % (int)(pow2_Nb_bits_minus_1);
        individus_bin[i].code_bin_deltaLambdaL = convert_to_binaire(deltaLambdaL);
        individus[i].deltaLambdaL = ((((deltaLambdaL_max - deltaLambdaL_min) / (pow2_Nb_bits_minus_1) * deltaLambdaL) + deltaLambdaL_min));
    }
}

/*________________________________________*/

long int Y0_en_bin(double Y0_dec)
{
    // Normalize the decimal value
    int Y0_ent;
    Y0_ent =  ((Y0_dec - Y0_min) * ((pow2_Nb_bits_minus_1)) / (Y0_max - Y0_min));
    
    // Convert the normalized value to binary
    return convert_to_binaire(Y0_ent);
}
/*________________________________________*/

long int B_en_bin(double B_dec)
{
    // Convert decimal value to integer in the range of [0, 2^Nb_bits - 1]
    int B_ent;
    B_ent = ((B_dec - B_min) * (pow2_Nb_bits_minus_1) / (B_max - B_min));

    // Convert integer to binary representation
    return convert_to_binaire(B_ent);
}
/*________________________________________*/

long int lambda0_en_bin(double lambda0_dec)
{
    // Calculate Lambda0 entity in 10-bit format
    int lambda0_ent;
    lambda0_ent = ((lambda0_dec - lambda0_min) * (pow2_Nb_bits_minus_1) / (lambda0_max - lambda0_min));
    
    // Convert to binary format
    return convert_to_binaire(lambda0_ent);
}
/*________________________________________*/

long int deltaLambdaL_en_bin(double deltaLambdaL_dec)
{
    // calculate and store integer deltaLambdaL_ent
    int deltaLambdaL_ent = ((deltaLambdaL_dec - deltaLambdaL_min) * (pow2_Nb_bits_minus_1) / (deltaLambdaL_max - deltaLambdaL_min));
    
    // convert deltaLambdaL_ent to binary format
    return convert_to_binaire(deltaLambdaL_ent);
}
/*________________________________________*/

void Ind_dec_to_Ind_bin(Individue_bin *ind_bin, Individue_dec *ind_dec)
{
    // Convert Y0 to binary
    ind_bin->code_bin_Y0 = Y0_en_bin(ind_dec->Y0);

    // Convert B to binary
    ind_bin->code_bin_B = B_en_bin(ind_dec->B);

    // Convert lambda0 to binary
    ind_bin->code_bin_lambda0 = lambda0_en_bin(ind_dec->lambda0);

    // Convert deltaLambdaL to binary
    ind_bin->code_bin_deltaLambdaL = deltaLambdaL_en_bin(ind_dec->deltaLambdaL);

    // Copy error code
    ind_bin->err = ind_dec->err;
}
/*________________________________________*/

// Generate array of unique integers
int *generer_nombres(int taille)
{
    // Allocate memory for array
    int *nombres = malloc(taille * sizeof(int));
    
    // Seed random number generator
    //srand(time(NULL));
    
    // Initialize array with unique values
    for (int i = 0; i < taille; i++)
    {
        nombres[i] = i;
    }

    // Shuffle array using Fisher-Yates algorithm
    for (int i = 0; i < taille; i++)
    {
        int j = rand() % taille;
        int temp = nombres[i];
        nombres[i] = nombres[j];
        nombres[j] = temp;
    }
    
    // Return shuffled array
    return nombres;
}
/*________________________________________*/

void Duel(Individue_bin individues[], Individue_bin nouvellePop[], int taille, int *new_pop_size)
{
    int *nombres = generer_nombres(taille); // generate indices
    int j = 0;
    for (int i = 0; i < taille; i += 2) // iterate over pairs of individuals
    {
        if (individues[nombres[i]].err < individues[nombres[i + 1]].err) // if first individual has smaller error
        {
            nouvellePop[j] = individues[nombres[i]]; // copy first individual to new population
        }
        else // if second individual has smaller error
        {
            nouvellePop[j] = individues[nombres[i + 1]]; // copy second individual to new population
        }
        j++; // increment index of new population
    }
    *new_pop_size = j;
    free(nombres); // free memory allocated for indices
}
/*________________________________________*/

int exposant(int base, int puissance)
{
    // Exponentiation function
    if (puissance == 0)
    {
        return 1;
    }
    else
    {
        int puiss = base;
        for (int i = 1; i < puissance; i++)
        {
            puiss *= base;
        }
        return puiss;
    }
}
/*________________________________________*/

void mutate_chromosome(long int *individu)
{
    int rando = rand() % 100; // Get random number

    // Check if mutation should occur
    if (rando < MUTATION_RATE)
    {
        // Select random bit to mutate
        int bit = rand() % (Nb_bits + 1);

        // Apply mutation
        int a = exposant(10, bit);
        if (((*individu / (a)) % 10) == 1)
        {
            *individu -= a; // Flip bit from 1 to 0
        }
        else
        {
            *individu += a; // Flip bit from 0 to 1
        }
    }
}
/*________________________________________*/

void mutate(Individue_bin *individu)
{
    // Mutate Y0 chromosome
    mutate_chromosome(&individu->code_bin_Y0);
    
    // Mutate B chromosome
    mutate_chromosome(&individu->code_bin_B);
    
    // Mutate lambda0 chromosome
    mutate_chromosome(&individu->code_bin_lambda0);
    
    // Mutate deltaLambdaL chromosome
    mutate_chromosome(&individu->code_bin_deltaLambdaL);
}
/*________________________________________*/

void croisement_chromosome(long int ind1, long int ind2, long int *enfant1, long int *enfant2)
{
    long int parent1_1ere, parent1_2eme, parent2_1ere, parent2_2eme;

    // Select a random bit for crossover
    int bit = rand() % (Nb_bits + 1);

    // Divide the parent chromosomes into two halves at the chosen bit
    parent1_1ere = (long int)ind1 % (long int)(exposant(10, (bit + 1)));
    parent1_2eme = (long int)ind1 / (long int)(exposant(10, (bit + 1)));
    parent2_1ere = (long int)ind2 % (long int)(exposant(10, (bit + 1)));
    parent2_2eme = (long int)ind2 / (long int)(exposant(10, (bit + 1)));

    // Swap the halves between the parent chromosomes to generate the children
    *enfant1 = parent1_1ere + (parent2_2eme * exposant(10, (bit + 1)));
    *enfant2 = parent2_1ere + (parent1_2eme * exposant(10, (bit + 1)));
}
/*________________________________________*/

void croisement(Individue_bin parent1, Individue_bin parent2, Individue_bin *enfant1, Individue_bin *enfant2)
{
    // Apply crossover on chromosomes
    croisement_chromosome(parent1.code_bin_Y0, parent2.code_bin_Y0, &enfant1->code_bin_Y0, &enfant2->code_bin_Y0);
    croisement_chromosome(parent1.code_bin_B, parent2.code_bin_B, &enfant1->code_bin_B, &enfant2->code_bin_B);
    croisement_chromosome(parent1.code_bin_lambda0, parent2.code_bin_lambda0, &enfant1->code_bin_lambda0, &enfant2->code_bin_lambda0);
    croisement_chromosome(parent1.code_bin_deltaLambdaL, parent2.code_bin_deltaLambdaL, &enfant1->code_bin_deltaLambdaL, &enfant2->code_bin_deltaLambdaL);
}
/*________________________________________*/

double densite(double deltaLambdaL)
{
    // Calcule de la densite
    double res;
    res = pow (10,21.78515 + 1.2614*log10(deltaLambdaL)) ;
    return res ;
}

/********************************************************************************************************************/