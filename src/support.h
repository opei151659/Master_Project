/*
	計算support set特徵值的過程
*/

#ifndef SUPPORT_H
#define SUPPORT_H

#include <iostream>
#include "tool.h"
extern "C" {
#include "espresso.h"   /* ESPRESSO.lib*/
}
#define low cube.last_word[cube.output]
#define fow cube.first_word[cube.output]
#define varnum BPI/2


pset_family support_set(pset_family F, pset_family R);
pcover sep_sup_output(pset_family T, int i, pset sup);

#endif SUPPORT_H