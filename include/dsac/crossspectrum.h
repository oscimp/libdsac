/*!
 * \file crossspectrum.h
 * \brief Description courte...

 * Description longue...
 */

#ifndef DSAC_CROSSSPECTRUM_H
#define DSAC_CROSSSPECTRUM_H

#ifdef __cplusplus /* if C++, specify external linkage to C functions */
extern "C" {
#endif

#include <stddef.h>

/*!
 * \brief Cross spectrum
 * Calcul le cross spectrum de deux voies complex.
 * La structure d'un chanel est un matrice à deux dimensions.
 * La première dimension sert à obtenir la partie réelle ou imaginaire d'un sample.
 * La seconde dimension contient tous les samples
 *
 * \param chanel_a Données du chanel A
 * \param chanel_b Données du chanel B
 * \param cross_spectrum Résultats du caclul
 * \param nb_elements Nombre de samples par chanel
 */
void cross_spectrum(const double * const chanel_a[], const double * const chanel_b[], double cross_spectrum[], const size_t nb_elements);

#ifdef __cplusplus
}
#endif

#endif // DSAC_CROSSSPECTRUM_H
