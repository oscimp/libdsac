#ifndef _FIR_H
#define _FIR_H

#include <stdlib.h>

#ifdef __cplusplus /* if C++, specify external linkage to C functions */
extern "C" {
#endif

/*!
 * \brief FIR version entier sur 32 bits et complexe
 * Calcul le FIR d'un jeu de données complexes codées sur 32 bits.
 *
 * \param vali Tableau contenant les i
 * \param valq Tableau contenant les y
 * \param data_size Taille des tableaux des données
 * \param coeff Tableau des coefficients
 * \param coeff_size Taille du tableau de coefficients
 * \param decim Facteur de decimation
 * \param accumi Resultat pour les i
 * \param accumq Resultat pour les q
 */
void fir_complex_32bits(const int32_t vali[], const int32_t valq[], const size_t data_size,
	const int32_t coeff[], const size_t coeff_size,
	const unsigned int decim, int64_t accumi[], int64_t accumq[]);

/*!
 * \brief FIR version double et complex
 * Calcul le FIR d'un jeu de données complexes codées en floatant double précision.
 *
 * \param vali Tableau contenant les i
 * \param valq Tableau contenant les y
 * \param data_size Taille des tableaux des données
 * \param coeff Tableau des coefficients
 * \param coeff_size Taille du tableau de coefficients
 * \param decim Facteur de decimation
 * \param accumi Resultat pour les i
 * \param accumq Resultat pour les q
 */
void fir_complex_double(const double vali[], const double valq[], const size_t data_size,
	const double coeff[], const size_t coeff_size, const unsigned long decim,
	double accumi[], double accumq[]);

/*!
 * \brief FIR version double et nombres réels
 * Calcul le FIR d'un jeu de données réel codées en floatant double précision.
 *
 * \param x Tableau contenant les données à filtrer
 * \param data_size Taille du tableau de données
 * \param coeff Tableau des coefficients
 * \param coeff_size Taille du tableau de coefficients
 * \param decim Facteur de decimation
 * \param accumx Contiendra les données filtrées
 */
void fir_double(const double x[], const size_t data_size, const double coeff[], const size_t coeff_size,
	const unsigned long decim, double accumx[]);

/*!
 * \brief FIR version float et nombres réels
 * Calcul le FIR d'un jeu de données réel codées en floatant simple précision.
 *
 * \param x Tableau contenant les données à filtrer
 * \param data_size Taille du tableau de données
 * \param coeff Tableau des coefficients
 * \param coeff_size Taille du tableau de coefficients
 * \param decim Facteur de decimation
 * \param accumx Contiendra les données filtrées
 */
void fir_float(const float x[], const size_t data_size, const float coeff[], const size_t coeff_size,
    const unsigned long decim, float accumx[]);

/*!
 * \brief FIR version entier signé 64 bits et nombres réels
 * Calcul le FIR d'un jeu de données réel codées en floatant double précision.
 *
 * \param x Tableau contenant les données à filtrer
 * \param data_size Taille du tableau de données
 * \param coeff Tableau des coefficients
 * \param coeff_size Taille du tableau de coefficients
 * \param decim Facteur de decimation
 * \param accumx Contiendra les données filtrées
 */
void fir_int64(const int64_t x[], const size_t data_size, const int64_t coeff[], const size_t coeff_size,
        const unsigned long decim,  int64_t accumx[]);

void fird(double *coeff, double *x,  int coeff_size, long data_size,
                int decim, double *accumx);
#ifdef __cplusplus
}
#endif

#endif // _FIR_H

